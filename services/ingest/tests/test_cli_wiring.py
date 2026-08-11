"""Guards on how the CLI stages call each other.

`run-all` chains the individual stage commands by calling them as plain Python
functions. That is convenient and has one sharp edge: a Typer command invoked
directly does NOT receive its declared defaults. Parameters left out arrive as
`typer.models.OptionInfo` objects, which sail through until something tries to
use the value.

That is not hypothetical -- `index()` was being called with no arguments, so
`run-all` died with "unknown source <OptionInfo object at 0x...>" *after* the
load had already completed. The documented one-command path was broken.

A plain "does run-all work" test would only catch it for the options that exist
today. This checks the property directly, so adding a new Option to any stage
fails here instead of in someone's multi-hour ingest.
"""

from __future__ import annotations

import ast
import inspect
import textwrap

import typer

from chemmed_ingest import cli

# Stage commands that run-all chains together.
STAGES = ["download", "parse", "fingerprint", "migrate", "load_cmd", "index"]


def option_parameters(fn) -> set[str]:
    """Parameters whose default is a Typer Option -- i.e. those that become
    OptionInfo objects when the function is called directly."""
    return {
        name
        for name, param in inspect.signature(fn).parameters.items()
        if isinstance(param.default, typer.models.OptionInfo)
    }


def calls_made_by(fn) -> dict[str, set[str]]:
    """Map of {called_function: {keyword argument names}} inside `fn`."""
    tree = ast.parse(textwrap.dedent(inspect.getsource(fn)))
    calls: dict[str, set[str]] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            calls.setdefault(node.func.id, set()).update(
                kw.arg for kw in node.keywords if kw.arg
            )
    return calls


class TestRunAllWiring:
    def test_run_all_passes_real_values_not_optioninfo(self):
        calls = calls_made_by(cli.run_all)

        problems: list[str] = []
        for stage in STAGES:
            if stage not in calls:
                continue
            required = option_parameters(getattr(cli, stage))
            missing = required - calls[stage]
            if missing:
                problems.append(f"{stage}() is missing {sorted(missing)}")

        assert not problems, (
            "run-all calls these stages without passing every option, so they "
            "will receive OptionInfo objects instead of values:\n  "
            + "\n  ".join(problems)
        )

    def test_run_all_actually_invokes_every_stage(self):
        """If a stage stops being called, the pipeline silently skips work."""
        calls = calls_made_by(cli.run_all)
        for stage in STAGES:
            assert stage in calls, f"run-all no longer calls {stage}()"

    def test_positional_args_are_not_used(self):
        """Positional calls would satisfy the check above while still being
        fragile to signature reordering."""
        tree = ast.parse(textwrap.dedent(inspect.getsource(cli.run_all)))
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
                and node.func.id in STAGES
            ):
                assert not node.args, f"{node.func.id}() called positionally"


class TestStageSignatures:
    def test_every_stage_exists(self):
        for stage in STAGES:
            assert callable(getattr(cli, stage, None)), f"{stage} is not defined"

    def test_index_source_option_is_still_named_source(self):
        """run-all passes source="parquet" by name; a rename would break it."""
        assert "source" in option_parameters(cli.index)

    def test_download_verify_option_is_still_named_verify(self):
        assert "verify" in option_parameters(cli.download)
