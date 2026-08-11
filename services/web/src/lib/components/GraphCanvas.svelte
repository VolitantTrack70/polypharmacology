<script lang="ts">
	import { onMount, onDestroy } from 'svelte';
	import type { GraphNode, GraphEdge } from '$lib/api';

	let { nodes = [], edges = [] }: { nodes: GraphNode[]; edges: GraphEdge[] } = $props();

	let container: HTMLDivElement;
	let cy: any = null;
	let selected = $state<GraphNode | null>(null);

	// Colour encodes node kind. The cascade reads left-to-right in the layout,
	// so the palette runs warm (query) to cool (pathway).
	const PALETTE: Record<string, string> = {
		query: '#e8590c',
		compound: '#f08c00',
		target: '#1971c2',
		pathway: '#2f9e44'
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
			],
			wheelSensitivity: 0.2
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
		background: #fbfbfc;
		border: 1px solid #e9ecef;
		border-radius: 8px;
		overflow: hidden;
	}
	.canvas {
		width: 100%;
		height: 100%;
	}
	.legend {
		position: absolute;
		bottom: 10px;
		left: 10px;
		display: flex;
		gap: 12px;
		font-size: 11px;
		color: #495057;
		background: rgba(255, 255, 255, 0.9);
		padding: 6px 10px;
		border-radius: 6px;
	}
	.chip {
		display: flex;
		align-items: center;
		gap: 5px;
	}
	.chip i {
		width: 9px;
		height: 9px;
		border-radius: 50%;
		display: inline-block;
	}
	.inspector {
		position: absolute;
		top: 10px;
		right: 10px;
		width: 240px;
		background: #fff;
		border: 1px solid #dee2e6;
		border-radius: 8px;
		padding: 10px 12px;
		font-size: 12px;
		box-shadow: 0 4px 16px rgba(0, 0, 0, 0.08);
	}
	.inspector header {
		display: flex;
		justify-content: space-between;
		align-items: start;
		gap: 8px;
		margin-bottom: 8px;
	}
	.inspector button {
		border: none;
		background: none;
		cursor: pointer;
		font-size: 16px;
		line-height: 1;
		color: #868e96;
	}
	dl {
		display: grid;
		grid-template-columns: auto 1fr;
		gap: 3px 10px;
		margin: 0;
	}
	dt {
		color: #868e96;
	}
	dd {
		margin: 0;
		word-break: break-word;
	}
	.empty {
		position: absolute;
		inset: 0;
		display: grid;
		place-content: center;
		color: #adb5bd;
		font-size: 13px;
		pointer-events: none;
	}
</style>
