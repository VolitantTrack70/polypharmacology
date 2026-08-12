<script lang="ts">
	import GraphCanvas from '$lib/components/GraphCanvas.svelte';
	import {
		api,
		ApiError,
		EXAMPLES,
		getDataSource,
		type DataSource,
		type OffTargetResponse
	} from '$lib/api';

	let smiles = $state('');
	let tanimoto = $state(0.4);
	let pchembl = $state(6.0);
	let organism = $state('Homo sapiens');

	let result = $state<OffTargetResponse | null>(null);
	let source = $state<DataSource>('unknown');
	let loading = $state(false);
	let error = $state<string | null>(null);
	let sortBy = $state<'affinity' | 'similarity'>('affinity');

	// Display filters. These never refetch -- the API returns the complete
	// result and the canvas shows a subset of it. The table stays complete.
	const QUERY_ID = '__query__';
	let topTargets = $state(25);
	let showCompounds = $state(false);
	let showPathways = $state(true);
	let targetFilter = $state('');

	// Index provenance in the header: what corpus a result was actually
	// computed against is part of the result.
	let indexInfo = $state<{ compounds: number | null; signature: string | null }>({
		compounds: null,
		signature: null
	});

	$effect(() => {
		api
			.status()
			.then((s) => {
				indexInfo = {
					compounds: (s.compounds_indexed as number) ?? null,
					signature: (s.fp_signature as string) ?? null
				};
			})
			.catch(() => {});
	});

	/** Best observed affinity per target, used to rank what to display. */
	let targetAffinity = $derived.by(() => {
		const m = new Map<string, number>();
		for (const e of result?.edges ?? []) {
			if (e.kind === 'binds_to') {
				m.set(e.target, Math.max(m.get(e.target) ?? 0, e.pchembl ?? 0));
			}
		}
		return m;
	});

	let rankedTargets = $derived.by(() => {
		const q = targetFilter.trim().toLowerCase();
		return (result?.nodes ?? [])
			.filter((n) => n.kind === 'target')
			.filter((n) => !q || n.label.toLowerCase().includes(q))
			.sort((a, b) => (targetAffinity.get(b.id) ?? 0) - (targetAffinity.get(a.id) ?? 0));
	});

	let displayGraph = $derived.by(() => {
		if (!result) return { nodes: [], edges: [] };

		const shown = topTargets > 0 ? rankedTargets.slice(0, topTargets) : rankedTargets;
		const keepTargets = new Set(shown.map((t) => t.id));

		const keepCompounds = new Set<string>();
		const keepPathways = new Set<string>();
		for (const e of result.edges) {
			if (showCompounds && e.kind === 'binds_to' && keepTargets.has(e.target)) {
				keepCompounds.add(e.source);
			}
			if (showPathways && e.kind === 'participates_in' && keepTargets.has(e.source)) {
				keepPathways.add(e.target);
			}
		}

		const keep = new Set([QUERY_ID, ...keepTargets, ...keepCompounds, ...keepPathways]);
		const nodes = result.nodes.filter((n) => keep.has(n.id));
		const edges = result.edges.filter((e) => keep.has(e.source) && keep.has(e.target));

		// With compounds hidden the query node would be orphaned, so collapse
		// query -> compound -> target into a direct edge carrying best affinity.
		if (!showCompounds) {
			for (const t of shown) {
				edges.push({
					source: QUERY_ID,
					target: t.id,
					kind: 'binds_to',
					pchembl: targetAffinity.get(t.id)
				});
			}
		}
		return { nodes, edges };
	});

	// 0.85 is the value intuition suggests and it returns essentially nothing
	// on ECFP4. Warn rather than block -- it's legitimate for near-duplicates.
	let thresholdWarning = $derived(
		tanimoto > 0.6
			? `At ${tanimoto.toFixed(2)}, ECFP4 matches only near-identical structures. Imatinib and nilotinib score 0.52.`
			: null
	);

	/** Flatten the node/edge graph into a target-centric table. */
	let rows = $derived.by(() => {
		if (!result) return [];
		const nodeById = new Map(result.nodes.map((n) => [n.id, n]));

		return result.nodes
			.filter((n) => n.kind === 'target')
			.map((target) => {
				const binding = result!.edges.filter(
					(e) => e.kind === 'binds_to' && e.target === target.id
				);
				const pathways = result!.edges
					.filter((e) => e.kind === 'participates_in' && e.source === target.id)
					.map((e) => nodeById.get(e.target)?.label)
					.filter((l): l is string => !!l);

				// Which similar compound gives the strongest evidence for this target.
				const via = binding
					.map((e) => ({
						compound: nodeById.get(e.source),
						pchembl: e.pchembl ?? 0,
						n: e.n_measurements ?? 0
					}))
					.sort((a, b) => b.pchembl - a.pchembl);

				return {
					id: target.id,
					gene: target.label,
					organism: target.organism,
					bestPchembl: Math.max(...binding.map((e) => e.pchembl ?? 0), 0),
					bestTanimoto: Math.max(...via.map((v) => v.compound?.tanimoto ?? 0), 0),
					viaCompound: via[0]?.compound?.label ?? '--',
					measurements: via.reduce((s, v) => s + v.n, 0),
					pathways
				};
			})
			.sort((a, b) =>
				sortBy === 'affinity'
					? b.bestPchembl - a.bestPchembl
					: b.bestTanimoto - a.bestTanimoto
			);
	});

	async function run() {
		if (!smiles.trim()) return;
		loading = true;
		error = null;
		try {
			result = await api.offTargets({
				smiles: smiles.trim(),
				tanimoto,
				pchembl,
				organism: organism || undefined
			});
		} catch (e) {
			error = e instanceof ApiError ? e.message : String(e);
			result = null;
		} finally {
			source = getDataSource();
			loading = false;
		}
	}

	function useExample(s: string) {
		smiles = s;
		run();
	}

	/** pChEMBL is -log10(molar). 6.0 -> 1000 nM. */
	function nanomolar(p: number): string {
		const nm = 10 ** (9 - p);
		if (nm >= 1000) return `${(nm / 1000).toFixed(1)} uM`;
		if (nm >= 1) return `${nm.toFixed(0)} nM`;
		return `${(nm * 1000).toFixed(0)} pM`;
	}
</script>

<div class="shell">
	<header class="topbar">
		<div>
			<h1>Polypharmacology &amp; Off-Target Graph</h1>
			<p class="mono sub">
				structure -> similar compounds -> known protein targets -> pathways
			</p>
		</div>
		<dl class="meta mono">
			<div>
				<dt>source</dt>
				<dd class:demo={source === 'demo'}>
					{source === 'demo' ? 'DEMO' : source === 'live' ? 'ChEMBL 35' : '--'}
				</dd>
			</div>
			<div>
				<dt>indexed</dt>
				<dd>{indexInfo.compounds ? indexInfo.compounds.toLocaleString() : '--'}</dd>
			</div>
			<div>
				<dt>fingerprint</dt>
				<dd>{indexInfo.signature ?? '--'}</dd>
			</div>
		</dl>
	</header>

	{#if source === 'demo'}
		<!--
			Non-negotiable. The fixtures are realistic enough to be mistaken for
			query output, and results that look real but aren't are worse than
			no results. Never remove this without removing the fallback.
		-->
		<div class="banner">
			<strong>Showing demo data.</strong> The API isn't running, so this is illustrative sample
			data -- real target biology, but the similarity and affinity numbers are not measurements.
			Start the backend to run genuine queries.
		</div>
	{/if}

	<div class="layout">
		<aside class="panel">
			<label class="field">
				<span class="lbl">Structure (SMILES)</span>
				<textarea
					bind:value={smiles}
					rows="3"
					spellcheck="false"
					placeholder="CC(=O)Oc1ccccc1C(=O)O"
					onkeydown={(e) => {
						if (e.key === 'Enter' && (e.metaKey || e.ctrlKey)) run();
					}}
				></textarea>
			</label>

			<div class="examples">
				{#each EXAMPLES as ex}
					<button class="chipbtn" onclick={() => useExample(ex.smiles)}>{ex.name}</button>
				{/each}
			</div>

			<label class="field">
				<span class="lbl">
					Tanimoto cutoff
					<b>{tanimoto.toFixed(2)}</b>
				</span>
				<input type="range" min="0.2" max="1" step="0.01" bind:value={tanimoto} />
				<div class="track-note">
					<i>0.20</i><i class="good">0.35-0.55</i><i>1.00</i>
				</div>
			</label>

			<label class="field">
				<span class="lbl">
					Min affinity (pChEMBL)
					<b>{pchembl.toFixed(1)}</b>
				</span>
				<input type="range" min="4" max="11" step="0.1" bind:value={pchembl} />
				<div class="track-note">
					<i>weaker</i><i class="good">{nanomolar(pchembl)}</i><i>potent</i>
				</div>
			</label>

			<label class="field">
				<span class="lbl">Organism</span>
				<select bind:value={organism}>
					<option value="Homo sapiens">Homo sapiens</option>
					<option value="">Any</option>
				</select>
			</label>

			{#if thresholdWarning}
				<p class="warn">{thresholdWarning}</p>
			{/if}

			<button class="primary" onclick={run} disabled={loading || !smiles.trim()}>
				{loading ? 'Searching...' : 'Find off-targets'}
			</button>
			<p class="kbd">Ctrl + Enter in the box also runs</p>
		</aside>

		<main class="main">
			{#if error}
				<p class="error">{error}</p>
			{/if}

			{#if result}
				<section class="stats">
					<div><b>{result.stats.similar_compounds}</b><span>similar compounds</span></div>
					<div><b>{result.stats.targets}</b><span>protein targets</span></div>
					<div><b>{result.stats.pathways}</b><span>pathways</span></div>
					<div>
						<b>{result.stats.search_ms.toFixed(0)}<small>ms</small></b>
						<span>{result.stats.compounds_scanned.toLocaleString()} scanned</span>
					</div>
				</section>

				{#if result.stats.truncated}
					<p class="hint">
						<strong>Showing a subset.</strong>
						{result.stats.similar_matched.toLocaleString()} compounds matched at Tanimoto &ge;
						{result.query.tanimoto_cutoff.toFixed(2)}, but the query limit cut it short -- some
						off-targets are missing. Raise the limit to see them all.
					</p>
				{/if}

				<section class="display">
					<div class="display-head">
						<h2>Graph display</h2>
						<span class="muted">
							{Math.min(topTargets || rankedTargets.length, rankedTargets.length)} of
							{targetFilter.trim()
								? `${rankedTargets.length} matching`
								: result.stats.targets} targets | {displayGraph.nodes.length} nodes
						</span>
					</div>
					<p class="muted note">
						Filters the canvas only. The table below always lists every target found.
					</p>

					<div class="display-controls">
						<label>
							<span>Top targets by affinity</span>
							<select bind:value={topTargets}>
								<option value={10}>10</option>
								<option value={25}>25</option>
								<option value={50}>50</option>
								<option value={0}>All ({result.stats.targets})</option>
							</select>
						</label>

						<label>
							<span>Filter by name</span>
							<input type="text" placeholder="e.g. kinase" bind:value={targetFilter} />
						</label>

						<label class="check">
							<input type="checkbox" bind:checked={showCompounds} />
							<span>Show similar compounds</span>
						</label>

						<label class="check">
							<input type="checkbox" bind:checked={showPathways} />
							<span>Show pathways</span>
						</label>
					</div>

					{#if rankedTargets.length === 0 && targetFilter.trim()}
						<p class="warn">No target matches "{targetFilter}".</p>
					{/if}
				</section>
			{/if}

			<section class="graph">
				<GraphCanvas nodes={displayGraph.nodes} edges={displayGraph.edges} />
			</section>

			{#if result && result.stats.similar_compounds === 0}
				<p class="hint">
					Nothing matched at Tanimoto &ge; {result.query.tanimoto_cutoff.toFixed(2)}. Try lowering
					it -- on ECFP4 fingerprints even closely related drugs typically score 0.4-0.5.
				</p>
			{/if}

			{#if rows.length}
				<section class="results">
					<div class="results-head">
						<h2>Predicted off-targets</h2>
						<div class="sort">
							<button class:on={sortBy === 'affinity'} onclick={() => (sortBy = 'affinity')}>
								by affinity
							</button>
							<button class:on={sortBy === 'similarity'} onclick={() => (sortBy = 'similarity')}>
								by similarity
							</button>
						</div>
					</div>

					<table>
						<thead>
							<tr>
								<th>Target</th>
								<th>Affinity</th>
								<th>Evidence via</th>
								<th>Tanimoto</th>
								<th>Pathways</th>
							</tr>
						</thead>
						<tbody>
							{#each rows as row (row.id)}
								<tr>
									<td>
										<b>{row.gene}</b>
										<small>{row.id}</small>
									</td>
									<td>
										<span class="mono">{row.bestPchembl.toFixed(1)}</span>
										<small>{nanomolar(row.bestPchembl)}</small>
									</td>
									<td>
										{row.viaCompound}
										<small>
											{row.measurements} measurement{row.measurements === 1 ? '' : 's'}
										</small>
									</td>
									<td>
										<div class="bar" style:--w="{row.bestTanimoto * 100}%">
											<span>{row.bestTanimoto.toFixed(2)}</span>
										</div>
									</td>
									<td class="pw">
										{#each row.pathways.slice(0, 2) as p}
											<span class="tag">{p}</span>
										{/each}
										{#if row.pathways.length > 2}
											<span class="more">+{row.pathways.length - 2}</span>
										{/if}
									</td>
								</tr>
							{/each}
						</tbody>
					</table>
				</section>
			{/if}
		</main>
	</div>

	<footer>
		<strong>These are hypotheses, not findings.</strong> Chemical similarity implies possible shared
		binding, not confirmed binding. Absence of a reported interaction is not evidence of absence --
		the underlying bioactivity data is heavily biased toward well-studied target families. Use this
		to decide what to test, never to conclude that something is safe.
	</footer>
</div>

<style>
	:global(:root) {
		--bg: #ffffff;
		--panel: #f4f4f2;
		--line: #b8bcc2;
		--line-soft: #dcdee1;
		--ink: #101317;
		--muted: #565c63;
		--accent: #14507a;
		--warn: #8a4b00;
		--good: #1f5c34;
		--danger: #8f1d1d;
		--mono: ui-monospace, 'Cascadia Mono', 'Consolas', 'DejaVu Sans Mono', monospace;
	}
	:global(body) {
		margin: 0;
		background: var(--bg);
		color: var(--ink);
		font-family: ui-sans-serif, system-ui, 'Segoe UI', sans-serif;
		font-size: 13px;
	}
	:global(*) {
		border-radius: 0;
	}
	.mono {
		font-family: var(--mono);
		font-variant-numeric: tabular-nums;
	}
	.shell {
		max-width: 1400px;
		margin: 0 auto;
		padding: 0 0 40px;
	}

	/* ---- header: title + index provenance ---- */
	.topbar {
		display: flex;
		justify-content: space-between;
		align-items: flex-end;
		gap: 24px;
		flex-wrap: wrap;
		padding: 14px 16px 10px;
		border-bottom: 2px solid var(--ink);
	}
	.topbar h1 {
		font-size: 15px;
		font-weight: 600;
		margin: 0;
		letter-spacing: -0.01em;
	}
	.sub {
		margin: 3px 0 0;
		font-size: 11px;
		color: var(--muted);
	}
	.meta {
		display: flex;
		gap: 0;
		margin: 0;
		border: 1px solid var(--line);
	}
	.meta > div {
		padding: 3px 10px;
		border-right: 1px solid var(--line);
	}
	.meta > div:last-child {
		border-right: 0;
	}
	.meta dt {
		font-size: 9.5px;
		text-transform: uppercase;
		letter-spacing: 0.06em;
		color: var(--muted);
	}
	.meta dd {
		margin: 1px 0 0;
		font-size: 11.5px;
	}
	.meta dd.demo {
		color: var(--warn);
		font-weight: 600;
	}

	.banner {
		margin: 0;
		padding: 7px 16px;
		background: #fdf6e8;
		border-bottom: 1px solid var(--line);
		border-left: 3px solid var(--warn);
		font-size: 12px;
		line-height: 1.5;
		color: #5c3200;
	}

	/* ---- layout ---- */
	.layout {
		display: grid;
		grid-template-columns: 260px 1fr;
		align-items: start;
	}
	@media (max-width: 900px) {
		.layout {
			grid-template-columns: 1fr;
		}
	}

	.panel {
		display: flex;
		flex-direction: column;
		gap: 12px;
		padding: 12px 14px;
		background: var(--panel);
		border-right: 1px solid var(--line);
		position: sticky;
		top: 0;
	}
	.field {
		display: flex;
		flex-direction: column;
		gap: 3px;
	}
	.lbl {
		font-size: 10px;
		text-transform: uppercase;
		letter-spacing: 0.06em;
		color: var(--muted);
		display: flex;
		justify-content: space-between;
	}
	.lbl b {
		font-family: var(--mono);
		font-variant-numeric: tabular-nums;
		color: var(--ink);
		text-transform: none;
		letter-spacing: 0;
	}
	textarea,
	select,
	input[type='text'] {
		font-family: var(--mono);
		font-size: 11.5px;
		padding: 5px 6px;
		border: 1px solid var(--line);
		background: #fff;
		color: var(--ink);
		resize: vertical;
	}
	textarea:focus,
	select:focus,
	input:focus {
		outline: 2px solid var(--accent);
		outline-offset: -2px;
	}
	input[type='range'] {
		width: 100%;
		accent-color: var(--accent);
		height: 16px;
	}
	.track-note {
		display: flex;
		justify-content: space-between;
		font-family: var(--mono);
		font-size: 9.5px;
		color: var(--muted);
	}
	.track-note i {
		font-style: normal;
	}
	.track-note .good {
		color: var(--good);
	}
	.examples {
		display: flex;
		flex-wrap: wrap;
		gap: 0;
		border: 1px solid var(--line);
	}
	.chipbtn {
		flex: 1 1 auto;
		font-size: 10.5px;
		padding: 4px 6px;
		border: 0;
		border-right: 1px solid var(--line);
		background: #fff;
		cursor: pointer;
		color: var(--muted);
	}
	.chipbtn:last-child {
		border-right: 0;
	}
	.chipbtn:hover {
		background: var(--accent);
		color: #fff;
	}
	.primary {
		padding: 7px;
		border: 1px solid var(--ink);
		background: var(--ink);
		color: #fff;
		font-size: 12px;
		font-weight: 500;
		cursor: pointer;
		letter-spacing: 0.02em;
	}
	.primary:hover:not(:disabled) {
		background: var(--accent);
		border-color: var(--accent);
	}
	.primary:disabled {
		background: #fff;
		color: #9aa0a6;
		border-color: var(--line);
		cursor: not-allowed;
	}
	.kbd {
		margin: -6px 0 0;
		font-size: 9.5px;
		color: var(--muted);
		text-align: center;
	}
	.warn {
		margin: 0;
		font-size: 11px;
		line-height: 1.45;
		color: var(--warn);
	}

	.main {
		display: flex;
		flex-direction: column;
		min-width: 0;
	}
	.error {
		margin: 0;
		padding: 8px 16px;
		background: #fdf0f0;
		border-bottom: 1px solid var(--danger);
		color: var(--danger);
		font-size: 12px;
	}

	/* ---- readout strip ---- */
	.stats {
		display: flex;
		flex-wrap: wrap;
		border-bottom: 1px solid var(--line);
	}
	.stats div {
		display: flex;
		flex-direction: column;
		padding: 8px 18px;
		border-right: 1px solid var(--line-soft);
	}
	.stats b {
		font-family: var(--mono);
		font-size: 19px;
		font-weight: 600;
		font-variant-numeric: tabular-nums;
		line-height: 1.1;
	}
	.stats small {
		font-size: 11px;
		font-weight: 400;
		color: var(--muted);
	}
	.stats span {
		font-size: 9.5px;
		text-transform: uppercase;
		letter-spacing: 0.06em;
		color: var(--muted);
		margin-top: 2px;
	}
	.hint {
		margin: 0;
		padding: 7px 16px;
		background: #fdf6e8;
		border-bottom: 1px solid var(--line);
		font-size: 12px;
		line-height: 1.5;
	}
	.canonical {
		display: block;
		font-family: var(--mono);
		font-size: 10.5px;
		color: var(--muted);
		padding: 5px 16px;
		border-bottom: 1px solid var(--line-soft);
		word-break: break-all;
	}

	/* ---- display filters ---- */
	.display {
		padding: 9px 16px;
		border-bottom: 1px solid var(--line);
		background: var(--panel);
	}
	.display-head {
		display: flex;
		justify-content: space-between;
		align-items: baseline;
		gap: 12px;
		flex-wrap: wrap;
	}
	.display h2 {
		font-size: 10px;
		font-weight: 600;
		text-transform: uppercase;
		letter-spacing: 0.07em;
		margin: 0;
		color: var(--muted);
	}
	.muted {
		color: var(--muted);
		font-size: 11px;
	}
	.display-head .muted {
		font-family: var(--mono);
		font-variant-numeric: tabular-nums;
	}
	.note {
		margin: 2px 0 8px;
		font-size: 10.5px;
	}
	.display-controls {
		display: flex;
		gap: 16px;
		flex-wrap: wrap;
		align-items: end;
	}
	.display-controls label {
		display: flex;
		flex-direction: column;
		gap: 3px;
		font-size: 10px;
		text-transform: uppercase;
		letter-spacing: 0.05em;
		color: var(--muted);
	}
	.display-controls select,
	.display-controls input[type='text'] {
		min-width: 140px;
	}
	.display-controls .check {
		flex-direction: row;
		align-items: center;
		gap: 5px;
		padding-bottom: 5px;
		cursor: pointer;
		text-transform: none;
		letter-spacing: 0;
		font-size: 11.5px;
	}
	.display-controls .check input {
		accent-color: var(--accent);
	}

	.graph {
		height: 520px;
		border-bottom: 1px solid var(--line);
	}

	/* ---- results table ---- */
	.results {
		padding: 10px 16px 0;
	}
	.results-head {
		display: flex;
		justify-content: space-between;
		align-items: baseline;
		margin-bottom: 6px;
	}
	.results h2 {
		font-size: 10px;
		font-weight: 600;
		text-transform: uppercase;
		letter-spacing: 0.07em;
		margin: 0;
		color: var(--muted);
	}
	.sort {
		display: flex;
		border: 1px solid var(--line);
	}
	.sort button {
		font-size: 10px;
		padding: 3px 9px;
		border: 0;
		border-right: 1px solid var(--line);
		background: #fff;
		cursor: pointer;
		color: var(--muted);
	}
	.sort button:last-child {
		border-right: 0;
	}
	.sort button.on {
		background: var(--ink);
		color: #fff;
	}
	table {
		width: 100%;
		border-collapse: collapse;
		font-size: 12px;
	}
	th {
		text-align: left;
		font-size: 9.5px;
		font-weight: 600;
		text-transform: uppercase;
		letter-spacing: 0.06em;
		color: var(--muted);
		padding: 5px 8px;
		border-top: 1px solid var(--ink);
		border-bottom: 1px solid var(--ink);
		background: var(--panel);
		white-space: nowrap;
	}
	td {
		padding: 5px 8px;
		border-bottom: 1px solid var(--line-soft);
		vertical-align: top;
	}
	tbody tr:nth-child(even) {
		background: #fafaf9;
	}
	tbody tr:hover {
		background: #eef3f7;
	}
	td b {
		font-weight: 600;
	}
	td small {
		display: block;
		font-family: var(--mono);
		font-size: 9.5px;
		color: var(--muted);
		margin-top: 1px;
	}
	.mono {
		font-family: var(--mono);
		font-variant-numeric: tabular-nums;
		font-weight: 600;
	}
	.bar {
		position: relative;
		background: #e6e8ea;
		border: 1px solid var(--line-soft);
		height: 14px;
		min-width: 58px;
	}
	.bar::before {
		content: '';
		position: absolute;
		inset: 0 auto 0 0;
		width: var(--w);
		background: #9fbcd2;
	}
	.bar span {
		position: relative;
		font-family: var(--mono);
		font-size: 10px;
		line-height: 14px;
		padding-left: 4px;
		font-variant-numeric: tabular-nums;
	}
	.pw {
		max-width: 300px;
	}
	.tag {
		display: inline-block;
		font-size: 9.5px;
		padding: 1px 5px;
		margin: 0 3px 2px 0;
		background: #eef2ef;
		border: 1px solid #c3d2c8;
		color: var(--good);
	}
	.more {
		font-family: var(--mono);
		font-size: 9.5px;
		color: var(--muted);
	}

	footer {
		margin: 22px 16px 0;
		padding-top: 10px;
		border-top: 1px solid var(--line);
		font-size: 10.5px;
		line-height: 1.6;
		color: var(--muted);
	}
</style>
