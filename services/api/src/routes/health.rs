use axum::extract::State;
use axum::Json;
use serde::Serialize;
use serde_json::{json, Value};

use crate::error::ApiResult;
use crate::AppState;

/// Liveness. Deliberately does not touch the database -- this answers "is the
/// process up", which is what a restart policy needs.
pub async fn health() -> Json<Value> {
    Json(json!({ "status": "ok", "service": "chemmed-api" }))
}

#[derive(Serialize)]
pub struct StatusResponse {
    database: ComponentStatus,
    chemworker: ComponentStatus,
    compounds_indexed: Option<u64>,
    fp_signature: Option<String>,
}

#[derive(Serialize)]
#[serde(tag = "state", rename_all = "snake_case")]
pub enum ComponentStatus {
    Up,
    Down { reason: String },
}

/// Readiness. Checks both dependencies and reports which one is broken --
/// "the graph is empty" and "the worker is down" look identical in the UI
/// otherwise.
pub async fn status(State(state): State<AppState>) -> ApiResult<Json<StatusResponse>> {
    let database = match sqlx::query_scalar::<_, i64>("SELECT 1")
        .fetch_one(&state.db)
        .await
    {
        Ok(_) => ComponentStatus::Up,
        Err(e) => ComponentStatus::Down { reason: e.to_string() },
    };

    let (chemworker, compounds_indexed, fp_signature) = match state.chem.status().await {
        Ok(s) => (ComponentStatus::Up, Some(s.compound_count), Some(s.fp_signature)),
        Err(e) => (ComponentStatus::Down { reason: e.message().to_string() }, None, None),
    };

    Ok(Json(StatusResponse { database, chemworker, compounds_indexed, fp_signature }))
}
