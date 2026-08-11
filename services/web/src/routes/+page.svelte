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
					viaCompound: via[0]?.compound?.label ?? '—',
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
		if (nm >= 1000) return `${(nm / 1000).toFixed(1)} µM`;
		if (nm >= 1) return `${nm.toFixed(0)} nM`;
		return `${(nm * 1000).toFixed(0)} pM`;
	}
</script>

<div class="shell">
	<header class="topbar">
		<div class="brand">
			<span class="mark"></span>
			<div>
				<h1>Polypharmacology &amp; Off-Target Graph</h1>
				<p>Structure &rarr; similar compounds &rarr; known protein targets &rarr; pathways</p>
			</div>
		</div>
		<span class="pill" class:demo={source === 'demo'} class:live={source === 'live'}>
			{source === 'demo' ? 'Demo data' : source === 'live' ? 'Live data' : 'Not connected'}
		</span>
	</header>

	{#if source === 'demo'}
		<!--
			Non-negotiable. The fixtures are realistic enough to be mistaken for
			query output, and results that look real but aren't are worse than
			no results. Never remove this without removing the fallback.
		-->
		<div class="banner">
			<strong>Showing demo data.</strong> The API isn't running, so this is illustrative sample
			data — real target biology, but the similarity and affinity numbers are not measurements.
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
					<i>0.20</i><i class="good">0.35–0.55</i><i>1.00</i>
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
				{loading ? 'Searching…' : 'Find off-targets'}
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
			{/if}

			<section class="graph">
				<GraphCanvas nodes={result?.nodes ?? []} edges={result?.edges ?? []} />
			</section>

			{#if result && result.stats.similar_compounds === 0}
				<p class="hint">
					Nothing matched at Tanimoto &ge; {result.query.tanimoto_cutoff.toFixed(2)}. Try lowering
					it — on ECFP4 fingerprints even closely related drugs typically score 0.4–0.5.
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
										<small>{row.measurements} measurements</small>
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
		binding, not confirmed binding. Absence of a reported interaction is not evidence of absence —
		the underlying bioactivity data is heavily biased toward well-studied target families. Use this
		to decide what to test, never to conclude that something is safe.
	</footer>
</div>

<style>
	:global(:root) {
		--bg: #ffffff;
		--panel: #fafbfc;
		--line: #e6e8eb;
		--ink: #14181d;
		--muted: #6b7480;
		--accent: #1f6feb;
		--warn: #d97706;
		--good: #16a34a;
	}
	:global(body) {
		margin: 0;
		background: var(--bg);
		color: var(--ink);
		font-family: ui-sans-serif, system-ui, -apple-system, 'Segoe UI', sans-serif;
		-webkit-font-smoothing: antialiased;
	}
	.shell {
		max-width: 1280px;
		margin: 0 auto;
		padding: 20px 22px 48px;
	}

	.topbar {
		display: flex;
		justify-content: space-between;
		align-items: center;
		gap: 16px;
		padding-bottom: 16px;
		border-bottom: 1px solid var(--line);
	}
	.brand {
		display: flex;
		gap: 12px;
		align-items: center;
	}
	.mark {
		width: 30px;
		height: 30px;
		border-radius: 8px;
		background: linear-gradient(135deg, #e8590c, #1f6feb 55%, #16a34a);
		flex: none;
	}
	.topbar h1 {
		font-size: 15px;
		margin: 0;
		letter-spacing: -0.01em;
	}
	.topbar p {
		margin: 2px 0 0;
		font-size: 11.5px;
		color: var(--muted);
	}
	.pill {
		font-size: 11px;
		padding: 4px 10px;
		border-radius: 99px;
		border: 1px solid var(--line);
		color: var(--muted);
		white-space: nowrap;
	}
	.pill.demo {
		background: #fff7ed;
		border-color: #fed7aa;
		color: #b45309;
	}
	.pill.live {
		background: #f0fdf4;
		border-color: #bbf7d0;
		color: #15803d;
	}

	.banner {
		margin-top: 14px;
		padding: 10px 14px;
		background: #fff7ed;
		border: 1px solid #fed7aa;
		border-left: 3px solid var(--warn);
		border-radius: 6px;
		font-size: 12.5px;
		line-height: 1.55;
		color: #7c2d12;
	}

	.layout {
		display: grid;
		grid-template-columns: 288px 1fr;
		gap: 22px;
		margin-top: 18px;
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
		gap: 15px;
		padding: 16px;
		background: var(--panel);
		border: 1px solid var(--line);
		border-radius: 10px;
		position: sticky;
		top: 16px;
	}
	.field {
		display: flex;
		flex-direction: column;
		gap: 6px;
	}
	.lbl {
		font-size: 11.5px;
		color: var(--muted);
		display: flex;
		justify-content: space-between;
	}
	.lbl b {
		color: var(--ink);
		font-variant-numeric: tabular-nums;
	}
	textarea,
	select {
		font: inherit;
		font-size: 11.5px;
		font-family: ui-monospace, 'Cascadia Code', Consolas, monospace;
		padding: 8px;
		border: 1px solid var(--line);
		border-radius: 6px;
		background: #fff;
		resize: vertical;
		color: var(--ink);
	}
	textarea:focus,
	select:focus {
		outline: 2px solid var(--accent);
		outline-offset: -1px;
	}
	input[type='range'] {
		width: 100%;
		accent-color: var(--accent);
	}
	.track-note {
		display: flex;
		justify-content: space-between;
		font-size: 10px;
		color: #9aa3ad;
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
		gap: 5px;
	}
	.chipbtn {
		font-size: 11px;
		padding: 3px 9px;
		border: 1px solid var(--line);
		background: #fff;
		border-radius: 99px;
		cursor: pointer;
		color: var(--muted);
	}
	.chipbtn:hover {
		border-color: var(--accent);
		color: var(--accent);
	}
	.primary {
		padding: 9px;
		border: none;
		border-radius: 6px;
		background: var(--accent);
		color: #fff;
		font-size: 12.5px;
		font-weight: 500;
		cursor: pointer;
	}
	.primary:disabled {
		background: #c3c9d0;
		cursor: not-allowed;
	}
	.kbd {
		margin: -8px 0 0;
		font-size: 10.5px;
		color: #9aa3ad;
		text-align: center;
	}
	.warn {
		margin: 0;
		font-size: 11.5px;
		line-height: 1.5;
		color: var(--warn);
	}

	.main {
		display: flex;
		flex-direction: column;
		gap: 14px;
		min-width: 0;
	}
	.error {
		margin: 0;
		padding: 10px 12px;
		background: #fef2f2;
		border: 1px solid #fecaca;
		border-radius: 6px;
		color: #b91c1c;
		font-size: 12.5px;
	}
	.stats {
		display: flex;
		gap: 30px;
		flex-wrap: wrap;
	}
	.stats div {
		display: flex;
		flex-direction: column;
	}
	.stats b {
		font-size: 21px;
		font-weight: 600;
		letter-spacing: -0.02em;
		font-variant-numeric: tabular-nums;
	}
	.stats small {
		font-size: 11px;
		font-weight: 400;
		color: var(--muted);
		margin-left: 1px;
	}
	.stats span {
		font-size: 10.5px;
		color: var(--muted);
	}
	.graph {
		height: 500px;
	}
	.hint {
		margin: 0;
		padding: 10px 12px;
		background: #fffbeb;
		border: 1px solid #fde68a;
		border-radius: 6px;
		font-size: 12.5px;
		line-height: 1.55;
	}

	.results-head {
		display: flex;
		justify-content: space-between;
		align-items: baseline;
		margin-bottom: 8px;
	}
	.results h2 {
		font-size: 13px;
		margin: 0;
	}
	.sort {
		display: flex;
		gap: 4px;
	}
	.sort button {
		font-size: 10.5px;
		padding: 3px 9px;
		border: 1px solid var(--line);
		background: #fff;
		border-radius: 99px;
		cursor: pointer;
		color: var(--muted);
	}
	.sort button.on {
		background: var(--ink);
		border-color: var(--ink);
		color: #fff;
	}
	table {
		width: 100%;
		border-collapse: collapse;
		font-size: 12px;
	}
	th {
		text-align: left;
		font-size: 10.5px;
		font-weight: 500;
		text-transform: uppercase;
		letter-spacing: 0.04em;
		color: var(--muted);
		padding: 7px 10px;
		border-bottom: 1px solid var(--line);
	}
	td {
		padding: 9px 10px;
		border-bottom: 1px solid #f1f3f5;
		vertical-align: top;
	}
	tbody tr:hover {
		background: #fafbfc;
	}
	td small {
		display: block;
		font-size: 10px;
		color: var(--muted);
		margin-top: 1px;
	}
	.mono {
		font-variant-numeric: tabular-nums;
		font-weight: 600;
	}
	.bar {
		position: relative;
		background: #eef1f4;
		border-radius: 3px;
		height: 16px;
		min-width: 54px;
	}
	.bar::before {
		content: '';
		position: absolute;
		inset: 0 auto 0 0;
		width: var(--w);
		background: #bfdbfe;
		border-radius: 3px;
	}
	.bar span {
		position: relative;
		font-size: 10.5px;
		line-height: 16px;
		padding-left: 5px;
		font-variant-numeric: tabular-nums;
	}
	.pw {
		max-width: 260px;
	}
	.tag {
		display: inline-block;
		font-size: 10px;
		padding: 2px 7px;
		margin: 0 3px 3px 0;
		background: #f0fdf4;
		border: 1px solid #bbf7d0;
		border-radius: 4px;
		color: #15803d;
	}
	.more {
		font-size: 10px;
		color: var(--muted);
	}

	footer {
		margin-top: 30px;
		padding-top: 14px;
		border-top: 1px solid var(--line);
		font-size: 11px;
		line-height: 1.65;
		color: var(--muted);
	}
</style>
