/**
 * Typed client for the Rust API.
 *
 * Shapes mirror `services/api/src/routes/offtargets.rs`. If you change one,
 * change both -- there is no codegen across this boundary yet.
 */

import { demoOffTargets, demoResolve } from './demo';

const BASE = import.meta.env.VITE_API_BASE ?? 'http://localhost:8080/api';

/**
 * When the real API is unreachable, fall back to demo fixtures so the UI is
 * always workable -- and always say so loudly.
 *
 * This exists because the backend needs Docker + Postgres + Rust + protoc to
 * come up, and the interface should not be blocked on any of that. `dataSource`
 * drives a persistent banner; results that look real but aren't are worse than
 * no results at all.
 */
export type DataSource = 'live' | 'demo' | 'unknown';

let dataSource: DataSource = 'unknown';
export const getDataSource = () => dataSource;

export type NodeKind = 'query' | 'compound' | 'target' | 'pathway';
export type EdgeKind = 'similar_to' | 'binds_to' | 'participates_in';

export interface GraphNode {
	id: string;
	label: string;
	kind: NodeKind;
	tanimoto?: number;
	organism?: string;
	biological_domain?: string;
}

export interface GraphEdge {
	source: string;
	target: string;
	kind: EdgeKind;
	tanimoto?: number;
	pchembl?: number;
	n_measurements?: number;
	activity_types?: string[];
}

export interface OffTargetResponse {
	query: {
		canonical_smiles: string;
		tanimoto_cutoff: number;
		pchembl_cutoff: number;
	};
	stats: {
		similar_compounds: number;
		targets: number;
		pathways: number;
		compounds_scanned: number;
		search_ms: number;
		/** Matches before `limit`. Greater than similar_compounds means a subset. */
		similar_matched: number;
		truncated: boolean;
	};
	nodes: GraphNode[];
	edges: GraphEdge[];
}

export interface ResolveResponse {
	canonical_smiles: string;
	inchikey: string;
	parent_inchikey: string;
	chembl_id: string | null;
	known_name: string | null;
}

export class ApiError extends Error {
	constructor(
		public readonly code: string,
		message: string
	) {
		super(message);
		this.name = 'ApiError';
	}
}

async function request<T>(path: string, init?: RequestInit): Promise<T> {
	let res: Response;
	try {
		res = await fetch(`${BASE}${path}`, {
			...init,
			headers: { 'Content-Type': 'application/json', ...init?.headers }
		});
	} catch {
		throw new ApiError('network', 'Could not reach the API. Is `cargo run` up?');
	}

	if (!res.ok) {
		// The API always returns { error: { code, message } }; fall back if it
		// somehow doesn't so the user sees something better than "[object Object]".
		const body = await res.json().catch(() => null);
		throw new ApiError(
			body?.error?.code ?? 'unknown',
			body?.error?.message ?? `Request failed (${res.status})`
		);
	}
	return res.json() as Promise<T>;
}

export interface OffTargetParams {
	smiles?: string;
	chembl_id?: string;
	tanimoto?: number;
	pchembl?: number;
	limit?: number;
	organism?: string;
	/** 'all' | 'domain' | 'specific' */
	pathway_scope?: string;
	max_pathways_per_target?: number;
	/** ChEMBL assay confidence floor, 7-9. */
	min_confidence?: number;
}

/**
 * Try the real API; fall back to fixtures if it isn't running.
 *
 * Only *connectivity* failures fall back. A 4xx from a live API is a genuine
 * result about the query (bad SMILES, no match) and must surface as an error --
 * silently swapping in demo data there would hide real bugs.
 */
async function withDemoFallback<T>(live: () => Promise<T>, demo: () => T): Promise<T> {
	try {
		const result = await live();
		dataSource = 'live';
		return result;
	} catch (e) {
		if (e instanceof ApiError && e.code === 'network') {
			dataSource = 'demo';
			return demo();
		}
		dataSource = 'live';
		throw e;
	}
}

export const api = {
	resolve: (smiles: string) =>
		withDemoFallback(
			() =>
				request<ResolveResponse>('/resolve', {
					method: 'POST',
					body: JSON.stringify({ smiles })
				}),
			() => demoResolve(smiles)
		),

	offTargets: (params: OffTargetParams) =>
		withDemoFallback(
			() =>
				request<OffTargetResponse>('/offtargets', {
					method: 'POST',
					body: JSON.stringify(params)
				}),
			() => demoOffTargets(params)
		),

	status: () => request<Record<string, unknown>>('/status')
};

/** A few structures worth having one click away while developing. */
export const EXAMPLES = [
	{ name: 'Imatinib', smiles: 'Cc1ccc(cc1Nc1nccc(n1)-c1cccnc1)NC(=O)c1ccc(cc1)CN1CCN(C)CC1' },
	{ name: 'Nilotinib', smiles: 'Cc1cn(cn1)-c1cc(cc(c1)C(F)(F)F)NC(=O)c1ccc(C)c(Nc2nccc(n2)-c2cccnc2)c1' },
	{ name: 'Aspirin', smiles: 'CC(=O)Oc1ccccc1C(=O)O' },
	{ name: 'Ibuprofen', smiles: 'CC(C)Cc1ccc(cc1)C(C)C(=O)O' }
];
