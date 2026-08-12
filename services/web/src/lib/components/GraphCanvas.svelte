<script lang="ts">
	import { onMount, onDestroy } from 'svelte';
	import type { GraphNode, GraphEdge } from '$lib/api';

	let { nodes = [], edges = [] }: { nodes: GraphNode[]; edges: GraphEdge[] } = $props();

	let container: HTMLDivElement;
	let cy: any = null;
	let selected = $state<GraphNode | null>(null);

	// Colour encodes node kind, desaturated to stay legible when printed
	// or projected rather than to look bright on screen.
	const PALETTE: Record<string, string> = {
		query: '#a03a00',
		compound: '#8a6d00',
		target: '#14507a',
		pathway: '#1f5c34'
	};

	onMount(async () => {
		const [{ default: cytoscape }, { default: fcose }] = await Promise.all([
			import('cytoscape'),
			import('cytoscape-fcose')
		]);
		cytoscape.use(fcose);

		cy = cytoscape({
			container,
			style: [
				{
					selector: 'node',
					style: {
						'background-color': (el: any) => PALETTE[el.data('kind')] ?? '#868e96',
						label: 'data(label)',
						'font-size': '9px',
						color: '#343a40',
						'text-valign': 'bottom',
						'text-margin-y': 4,
						// Long protein and pathway names otherwise swamp the canvas.
						'text-max-width': '110px',
						'text-wrap': 'ellipsis',
						width: (el: any) => (el.data('kind') === 'query' ? 34 : 18),
						height: (el: any) => (el.data('kind') === 'query' ? 34 : 18)
					}
				},
				{
					selector: 'edge',
					style: {
						'curve-style': 'bezier',
						'line-color': '#ced4da',
						width: (el: any) =>
							// Similarity edges carry a score worth seeing; affinity
							// edges encode pChEMBL the same way.
							el.data('kind') === 'similar_to'
								? Math.max(1, (el.data('tanimoto') ?? 0.4) * 5)
								: Math.max(1, ((el.data('pchembl') ?? 6) - 5) * 1.2),
						'target-arrow-shape': 'triangle',
						'target-arrow-color': '#ced4da',
						'arrow-scale': 0.6
					}
				},
				{
					selector: 'node:selected',
					style: { 'border-width': 3, 'border-color': '#212529' }
				}
			]
			// wheelSensitivity is deliberately left at the default -- overriding it
			// makes zoom behave unpredictably across mice and cytoscape warns.
		});

		cy.on('tap', 'node', (evt: any) => {
			selected = evt.target.data();
		});
		cy.on('tap', (evt: any) => {
			if (evt.target === cy) selected = null;
		});

		render();
	});

	onDestroy(() => cy?.destroy());

	function render() {
		if (!cy) return;
		cy.elements().remove();
		cy.add([
			...nodes.map((n) => ({ group: 'nodes' as const, data: { ...n } })),
			...edges.map((e) => ({
				group: 'edges' as const,
				data: { ...e, id: `${e.source}->${e.target}:${e.kind}` }
			}))
		]);
		if (nodes.length === 0) return;

		cy.layout({
			name: 'fcose',
			animate: false,
			quality: 'default',
			nodeRepulsion: 8000,
			idealEdgeLength: 70,
			randomize: true
		}).run();
		cy.fit(undefined, 40);
	}

	// Re-render whenever the parent hands us a new result set.
	$effect(() => {
		void nodes;
		void edges;
		render();
	});
</script>

<div class="wrap">
	<div class="canvas" bind:this={container}></div>

	<div class="legend">
		{#each Object.entries(PALETTE) as [kind, color]}
			<span class="chip"><i style:background={color}></i>{kind}</span>
		{/each}
	</div>

	{#if selected}
		<aside class="inspector">
			<header>
				<strong>{selected.label}</strong>
				<button onclick={() => (selected = null)} aria-label="Close">&times;</button>
			</header>
			<dl>
				<dt>Kind</dt>
				<dd>{selected.kind}</dd>
				<dt>ID</dt>
				<dd><code>{selected.id}</code></dd>
				{#if selected.tanimoto != null}
					<dt>Tanimoto</dt>
					<dd>{selected.tanimoto.toFixed(3)}</dd>
				{/if}
				{#if selected.organism}
					<dt>Organism</dt>
					<dd>{selected.organism}</dd>
				{/if}
				{#if selected.biological_domain}
					<dt>Domain</dt>
					<dd>{selected.biological_domain}</dd>
				{/if}
			</dl>
		</aside>
	{/if}

	{#if nodes.length === 0}
		<p class="empty">No cascade to show yet. Enter a structure and run a query.</p>
	{/if}
</div>

<style>
	.wrap {
		position: relative;
		width: 100%;
		height: 100%;
		min-height: 480px;
		background: #fdfdfc;
		overflow: hidden;
	}
	.canvas {
		width: 100%;
		height: 100%;
	}
	.legend {
		position: absolute;
		bottom: 0;
		left: 0;
		display: flex;
		font-size: 10px;
		color: #565c63;
		background: #fff;
		border-top: 1px solid #b8bcc2;
		border-right: 1px solid #b8bcc2;
	}
	.chip {
		display: flex;
		align-items: center;
		gap: 5px;
		padding: 3px 9px;
		border-right: 1px solid #dcdee1;
		text-transform: uppercase;
		letter-spacing: 0.05em;
	}
	.chip:last-child {
		border-right: 0;
	}
	.chip i {
		width: 8px;
		height: 8px;
		display: inline-block;
	}
	.inspector {
		position: absolute;
		top: 0;
		right: 0;
		width: 250px;
		background: #fff;
		border-left: 1px solid #b8bcc2;
		border-bottom: 1px solid #b8bcc2;
		padding: 8px 10px;
		font-size: 11.5px;
	}
	.inspector header {
		display: flex;
		justify-content: space-between;
		align-items: start;
		gap: 8px;
		margin-bottom: 6px;
		padding-bottom: 5px;
		border-bottom: 1px solid #dcdee1;
	}
	.inspector button {
		border: none;
		background: none;
		cursor: pointer;
		font-size: 15px;
		line-height: 1;
		color: #565c63;
	}
	dl {
		display: grid;
		grid-template-columns: auto 1fr;
		gap: 2px 10px;
		margin: 0;
	}
	dt {
		color: #565c63;
		font-size: 9.5px;
		text-transform: uppercase;
		letter-spacing: 0.05em;
		padding-top: 1px;
	}
	dd {
		margin: 0;
		word-break: break-word;
		font-family: ui-monospace, 'Cascadia Mono', Consolas, monospace;
		font-size: 10.5px;
	}
	.empty {
		position: absolute;
		inset: 0;
		display: grid;
		place-content: center;
		color: #8b9096;
		font-size: 12px;
		pointer-events: none;
	}
</style>
