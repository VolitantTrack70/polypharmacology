use axum::extract::{Path, State};
use axum::Json;
use serde::{Deserialize, Serialize};

use crate::error::{ApiError, ApiResult};
use crate::AppState;

#[derive(Deserialize)]
pub struct ResolveRequest {
    pub smiles: String,
}

#[derive(Serialize)]
pub struct ResolveResponse {
    pub canonical_smiles: String,
    pub inchikey: String,
    pub parent_inchikey: String,
    /// Populated when this exact structure is already in ChEMBL. Lets the UI
    /// say "this is CHEMBL941 (imatinib)" instead of treating a known drug as
    /// an anonymous query structure.
    pub chembl_id: Option<String>,
    pub known_name: Option<String>,
}

/// Standardise a user-supplied SMILES and look for an exact ChEMBL match.
///
/// Matching is on InChIKey rather than SMILES string: the same molecule has
/// many valid SMILES, and only the InChIKey is canonical across toolkits.
pub async fn resolve(
    State(state): State<AppState>,
    Json(req): Json<ResolveRequest>,
) -> ApiResult<Json<ResolveResponse>> {
    if req.smiles.trim().is_empty() {
        return Err(ApiError::BadRequest("smiles must not be empty".into()));
    }

    let std = state.chem.standardize(&req.smiles).await?;
    if !std.ok {
        return Err(ApiError::InvalidStructure(std.error));
    }

    let hit: Option<(String, Option<String>)> = sqlx::query_as(
        "SELECT chembl_id, pref_name FROM chem.compound \
         WHERE standard_inchi_key = $1 LIMIT 1",
    )
    .bind(&std.inchikey)
    .fetch_optional(&state.db)
    .await?;

    let (chembl_id, known_name) = match hit {
        Some((id, name)) => (Some(id), name),
        None => (None, None),
    };

    Ok(Json(ResolveResponse {
        canonical_smiles: std.canonical_smiles,
        inchikey: std.inchikey,
        parent_inchikey: std.parent_inchikey,
        chembl_id,
        known_name,
    }))
}

#[derive(Serialize, sqlx::FromRow)]
pub struct CompoundDetail {
    pub chembl_id: String,
    pub canonical_smiles: Option<String>,
    pub standard_inchi_key: Option<String>,
    pub molformula: Option<String>,
    pub mw_freebase: Option<f32>,
    pub alogp: Option<f32>,
    pub hba: Option<i16>,
    pub hbd: Option<i16>,
    pub psa: Option<f32>,
    pub num_ro5_violations: Option<i16>,
    pub max_phase: Option<f32>,
    pub first_approval: Option<i16>,
    pub withdrawn_flag: bool,
}

pub async fn get_compound(
    State(state): State<AppState>,
    Path(chembl_id): Path<String>,
) -> ApiResult<Json<CompoundDetail>> {
    let detail: Option<CompoundDetail> = sqlx::query_as(
        r#"
        SELECT chembl_id, canonical_smiles, standard_inchi_key, molformula,
               mw_freebase, alogp, hba, hbd, psa, num_ro5_violations,
               max_phase, first_approval, withdrawn_flag
        FROM chem.compound
        WHERE chembl_id = $1
        "#,
    )
    .bind(&chembl_id)
    .fetch_optional(&state.db)
    .await?;

    detail
        .map(Json)
        .ok_or_else(|| ApiError::NotFound(format!("no compound {chembl_id}")))
}
