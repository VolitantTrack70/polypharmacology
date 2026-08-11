/**
 * DEMO FIXTURES — synthetic data for developing the UI without a backend.
 *
 * ============================ READ THIS ============================
 * This is NOT real query output. It is hand-built sample data so the
 * interface can be designed, reviewed, and demonstrated before the
 * ingestion pipeline and API are running.
 *
 * The biology is real and well-documented (imatinib's kinase target
 * profile is textbook, the UniProt accessions are genuine, the Reactome
 * pathway names are real). The *numbers* — Tanimoto coefficients and
 * pChEMBL values — are plausible illustrative values, NOT measurements.
 *
 * Whenever this module supplies the data, the UI shows a persistent
 * "Demo data" banner. Never remove that banner: results that look real
 * but aren't are worse than no results.
 * ===================================================================
 */

import type { GraphEdge, GraphNode, OffTargetResponse } from './api';

/** Imatinib's documented target profile. Genuine gene symbols + accessions. */
const IMATINIB_TARGETS = [
	{ id: 'CHEMBL1862', gene: 'ABL1', uniprot: 'P00519', name: 'Tyrosine-protein kinase ABL1' },
	{ id: 'CHEMBL1936', gene: 'KIT', uniprot: 'P10721', name: 'Stem cell growth factor receptor' },
	{ id: 'CHEMBL2007', gene: 'PDGFRA', uniprot: 'P16234', name: 'PDGF receptor alpha' },
	{ id: 'CHEMBL1913', gene: 'PDGFRB', uniprot: 'P09619', name: 'PDGF receptor beta' },
	{ id: 'CHEMBL5122', gene: 'DDR1', uniprot: 'Q08345', name: 'Epithelial discoidin domain receptor 1' },
	{ id: 'CHEMBL3145', gene: 'NQO2', uniprot: 'P16083', name: 'Ribosyldihydronicotinamide dehydrogenase' }
];

/** Real Reactome pathway names. */
const PATHWAYS = [
	{ id: 'R-HSA-1433557', name: 'Signaling by SCF-KIT', domain: 'Signal Transduction' },
	{ id: 'R-HSA-186797', name: 'Signaling by PDGF', domain: 'Signal Transduction' },
	{ id: 'R-HSA-9006934', name: 'Signaling by Receptor Tyrosine Kinases', domain: 'Signal Transduction' },
	{ id: 'R-HSA-1236394', name: 'Signaling by ERBB2', domain: 'Signal Transduction' },
	{ id: 'R-HSA-8878171', name: 'Transcriptional regulation by RUNX1', domain: 'Gene expression' },
	{ id: 'R-HSA-9020702', name: 'Interleukin-1 signaling', domain: 'Immune System' }
];

/** Structurally related BCR-ABL inhibitors, with illustrative similarity. */
const ANALOGUES = [
	{ id: 'CHEMBL941', name: 'Imatinib', tanimoto: 1.0 },
	{ id: 'CHEMBL255863', name: 'Nilotinib', tanimoto: 0.517 },
	{ id: 'CHEMBL1421', name: 'Dasatinib', tanimoto: 0.431 },
	{ id: 'CHEMBL288441', name: 'Bosutinib', tanimoto: 0.402 },
	{ id: 'CHEMBL1171837', name: 'Ponatinib', tanimoto: 0.388 }
];

const TARGETS_PER_ANALOGUE: Record<string, string[]> = {
	CHEMBL941: ['ABL1', 'KIT', 'PDGFRA', 'PDGFRB', 'DDR1', 'NQO2'],
	CHEMBL255863: ['ABL1', 'KIT', 'PDGFRA', 'PDGFRB', 'DDR1'],
	CHEMBL1421: ['ABL1', 'KIT', 'PDGFRB'],
	CHEMBL288441: ['ABL1', 'KIT'],
	CHEMBL1171837: ['ABL1', 'PDGFRA', 'DDR1']
};

const PATHWAYS_PER_TARGET: Record<string, string[]> = {
	ABL1: ['R-HSA-9006934', 'R-HSA-8878171'],
	KIT: ['R-HSA-1433557', 'R-HSA-9006934'],
	PDGFRA: ['R-HSA-186797', 'R-HSA-9006934'],
	PDGFRB: ['R-HSA-186797', 'R-HSA-9006934'],
	DDR1: ['R-HSA-9006934'],
	NQO2: ['R-HSA-9020702']
};

/** Illustrative affinities (pChEMBL-like). Not measurements. */
const AFFINITY: Record<string, number> = {
	ABL1: 8.4, KIT: 7.9, PDGFRA: 7.6, PDGFRB: 7.8, DDR1: 7.1, NQO2: 6.3
};

const QUERY_ID = '__query__';

/**
 * Build a demo cascade. Honours the threshold sliders so the controls behave
 * correctly during UI work — dragging the cutoff genuinely changes the result.
 */
export function demoOffTargets(params: {
	smiles?: string;
	tanimoto?: number;
	pchembl?: number;
}): OffTargetResponse {
	const tCut = params.tanimoto ?? 0.4;
	const pCut = params.pchembl ?? 6.0;

	const analogues = ANALOGUES.filter((a) => a.tanimoto >= tCut);

	const nodes: GraphNode[] = [
		{ id: QUERY_ID, label: 'query structure', kind: 'query' }
	];
	const edges: GraphEdge[] = [];

	const seenTargets = new Set<string>();
	const seenPathways = new Set<string>();

	for (const analogue of analogues) {
		nodes.push({
			id: analogue.id,
			label: analogue.name,
			kind: 'compound',
			tanimoto: analogue.tanimoto
		});
		edges.push({
			source: QUERY_ID,
			target: analogue.id,
			kind: 'similar_to',
			tanimoto: analogue.tanimoto
		});

		for (const gene of TARGETS_PER_ANALOGUE[analogue.id] ?? []) {
			const affinity = AFFINITY[gene] ?? 6.5;
			if (affinity < pCut) continue;

			const target = IMATINIB_TARGETS.find((t) => t.gene === gene);
			if (!target) continue;

			if (!seenTargets.has(target.id)) {
				seenTargets.add(target.id);
				nodes.push({
					id: target.id,
					label: target.gene,
					kind: 'target',
					organism: 'Homo sapiens'
				});
			}

			edges.push({
				source: analogue.id,
				target: target.id,
				kind: 'binds_to',
				pchembl: affinity,
				n_measurements: Math.round(affinity * 7),
				activity_types: ['IC50', 'Kd']
			});

			for (const pid of PATHWAYS_PER_TARGET[gene] ?? []) {
				const pathway = PATHWAYS.find((p) => p.id === pid);
				if (!pathway) continue;

				if (!seenPathways.has(pathway.id)) {
					seenPathways.add(pathway.id);
					nodes.push({
						id: pathway.id,
						label: pathway.name,
						kind: 'pathway',
						biological_domain: pathway.domain
					});
				}
				if (!edges.some((e) => e.source === target.id && e.target === pathway.id)) {
					edges.push({
						source: target.id,
						target: pathway.id,
						kind: 'participates_in'
					});
				}
			}
		}
	}

	return {
		query: {
			canonical_smiles:
				params.smiles ?? 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1',
			tanimoto_cutoff: tCut,
			pchembl_cutoff: pCut
		},
		stats: {
			similar_compounds: analogues.length,
			targets: seenTargets.size,
			pathways: seenPathways.size,
			compounds_scanned: 2_409_270,
			search_ms: 143
		},
		nodes,
		edges
	};
}

export function demoResolve(smiles: string) {
	return {
		canonical_smiles: smiles,
		inchikey: 'KTUFNOKKBVMGRW-UHFFFAOYSA-N',
		parent_inchikey: 'KTUFNOKKBVMGRW',
		chembl_id: 'CHEMBL941',
		known_name: 'Imatinib'
	};
}
