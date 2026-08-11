<script lang="ts">
	import GraphCanvas from '$lib/components/GraphCanvas.svelte';
	import { api, ApiError, EXAMPLES, type OffTargetResponse } from '$lib/api';

	let smiles = $state('');
	let tanimoto = $state(0.4);
	let pchembl = $state(6.0);
	let organism = $state('Homo sapiens');

	let result = $state<OffTargetResponse | null>(null);
	let loading = $state(false);
	let error = $state<string | null>(null);

	// The 0.85 that intuition suggests returns essentially nothing on ECFP4.
	// Warn rather than block -- it is a legitimate setting for near-duplicates.
	let thresholdWarning = $derived(
		tanimoto > 0.6
			? 'Above ~0.6, ECFP4 matches only near-identical structures. Imatinib/nilotinib score 0.52.'
			: null
	);

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
			loading = false;
		}
	}

	function useExample(s: string) {
		smiles = s;
		run();
	}
</script>

<main>
	<header class="masthead">
		<h1>Polypharmacology &amp; Off-Target Graph</h1>
		<p class="sub">
			Finds structurally similar compounds, the proteins they are known to bind, and the
			pathways those proteins sit in.
		</p>
	</header>

	<section class="controls">
		<label class="field">
			<span>Structure (SMILES)</span>
			<textarea
				bind:value={smiles}
				rows="2"
				placeholder="CC(=O)Oc1ccccc1C(=O)O"
				onkeydown={(e) => {
					if (e.key === 'Enter' && (e.metaKey || e.ctrlKey)) run();
				}}
			></textarea>
		</label>

		<div class="examples">
			<span>Try:</span>
			{#each EXAMPLES as ex}
				<button class="link" onclick={() => useExample(ex.smiles)}>{ex.name}</button>
			{/each}
		</div>

		<div class="sliders">
			<label class="field">
				<span>
					Tanimoto cutoff <strong>{tanimoto.toFixed(2)}</strong>
				</span>
				<input type="range" min="0.2" max="1" step="0.01" bind:value={tanimoto} />
				<small class="scale"><i>0.20</i><i class="useful">0.35–0.55 useful</i><i>1.00</i></small>
			</label>

			<label class="field">
				<span>Min pChEMBL <strong>{pchembl.toFixed(1)}</strong></span>
				<input type="range" min="4" max="11" step="0.1" bind:value={pchembl} />
				<small class="scale"><i>weak</i><i>{(10 ** (9 - pchembl)).toFixed(0)} nM</i><i>potent</i></small>
			</label>

			<label class="field">
				<span>Organism</span>
				<select bind:value={organism}>
					<option value="Homo sapiens">Homo sapiens</option>
					<option value="">Any</option>
				</select>
			</label>
		</div>

		{#if thresholdWarning}
			<p class="warn">{thresholdWarning}</p>
		{/if}

		<button class="primary" onclick={run} disabled={loading || !smiles.trim()}>
			{loading ? 'Searching…' : 'Find off-targets'}
		</button>
	</section>

	{#if error}
		<p class="error">{error}</p>
	{/if}

	{#if result}
		<section class="stats">
			<div><strong>{result.stats.similar_compounds}</strong><span>similar compounds</span></div>
			<div><strong>{result.stats.targets}</strong><span>targets</span></div>
			<div><strong>{result.stats.pathways}</strong><span>pathways</span></div>
			<div>
				<strong>{result.stats.search_ms.toFixed(0)} ms</strong>
				<span>over {result.stats.compounds_scanned.toLocaleString()} compounds</span>
			</div>
		</section>

		{#if result.stats.similar_compounds === 0}
			<p class="hint">
				Nothing matched at Tanimoto ≥ {result.query.tanimoto_cutoff.toFixed(2)}. Try lowering the
				cutoff — on ECFP4 fingerprints, even closely related drugs often score around 0.4–0.5.
			</p>
		{/if}

		<code class="canonical">{result.query.canonical_smiles}</code>
	{/if}

	<section class="graph">
		<GraphCanvas nodes={result?.nodes ?? []} edges={result?.edges ?? []} />
	</section>

	<footer>
		<p>
			<strong>These are hypotheses, not findings.</strong> Chemical similarity implies possible
			shared binding, not confirmed binding. Absence of a reported interaction is not evidence of
			absence — the underlying bioactivity data is heavily biased toward well-studied target
			families. Use this to decide what to test, never to conclude that something is safe.
		</p>
	</footer>
</main>

<style>
	:global(body) {
		margin: 0;
		font-family: ui-sans-serif, system-ui, -apple-system, 'Segoe UI', sans-serif;
		color: #212529;
		background: #fff;
	}
	main {
		max-width: 1100px;
		margin: 0 auto;
		padding: 28px 20px 60px;
	}
	.masthead h1 {
		font-size: 22px;
		margin: 0 0 4px;
	}
	.sub {
		margin: 0 0 24px;
		color: #868e96;
		font-size: 13px;
	}
	.controls {
		display: flex;
		flex-direction: column;
		gap: 14px;
		padding: 18px;
		border: 1px solid #e9ecef;
		border-radius: 10px;
		background: #fbfbfc;
	}
	.field {
		display: flex;
		flex-direction: column;
		gap: 5px;
		font-size: 12px;
		color: #495057;
	}
	textarea,
	select {
		font: inherit;
		font-family: ui-monospace, 'Cascadia Code', monospace;
		font-size: 12px;
		padding: 8px;
		border: 1px solid #dee2e6;
		border-radius: 6px;
		resize: vertical;
	}
	.sliders {
		display: grid;
		grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
		gap: 16px;
	}
	.scale {
		display: flex;
		justify-content: space-between;
		font-size: 10px;
		color: #adb5bd;
		font-style: normal;
	}
	.scale i {
		font-style: normal;
	}
	.scale .useful {
		color: #2f9e44;
	}
	.examples {
		display: flex;
		gap: 8px;
		align-items: center;
		font-size: 12px;
		color: #868e96;
		flex-wrap: wrap;
	}
	.link {
		background: none;
		border: none;
		padding: 0;
		color: #1971c2;
		cursor: pointer;
		font-size: 12px;
		text-decoration: underline;
	}
	.primary {
		align-self: start;
		padding: 9px 18px;
		border: none;
		border-radius: 6px;
		background: #1971c2;
		color: #fff;
		font-size: 13px;
		cursor: pointer;
	}
	.primary:disabled {
		background: #adb5bd;
		cursor: not-allowed;
	}
	.warn {
		margin: 0;
		font-size: 12px;
		color: #e8590c;
	}
	.error {
		margin: 16px 0 0;
		padding: 10px 12px;
		border-radius: 6px;
		background: #fff5f5;
		color: #c92a2a;
		font-size: 13px;
	}
	.stats {
		display: flex;
		gap: 28px;
		margin: 22px 0 10px;
		flex-wrap: wrap;
	}
	.stats div {
		display: flex;
		flex-direction: column;
	}
	.stats strong {
		font-size: 20px;
	}
	.stats span {
		font-size: 11px;
		color: #868e96;
	}
	.hint {
		font-size: 13px;
		color: #495057;
		background: #fff9db;
		padding: 10px 12px;
		border-radius: 6px;
	}
	.canonical {
		display: block;
		font-size: 11px;
		color: #868e96;
		margin-bottom: 12px;
		word-break: break-all;
	}
	.graph {
		height: 560px;
		margin-top: 10px;
	}
	footer p {
		margin-top: 28px;
		font-size: 11px;
		line-height: 1.6;
		color: #868e96;
		border-top: 1px solid #e9ecef;
		padding-top: 14px;
	}
</style>
