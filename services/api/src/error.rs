//! Error type shared across handlers.
//!
//! Internal failures are logged with full detail but reported to the client as
//! a generic message -- database errors leak schema, and this API is meant to
//! be exposable beyond localhost eventually.

use axum::http::StatusCode;
use axum::response::{IntoResponse, Response};
use axum::Json;
use serde_json::json;

#[derive(Debug, thiserror::Error)]
pub enum ApiError {
    #[error("{0}")]
    BadRequest(String),

    #[error("{0}")]
    NotFound(String),

    /// The structure could not be parsed or standardised. Distinct from a
    /// generic bad request because the UI surfaces it inline on the input.
    #[error("invalid structure: {0}")]
    InvalidStructure(String),

    #[error("chemworker unavailable: {0}")]
    ChemWorker(#[from] tonic::Status),

    #[error("database error: {0}")]
    Database(#[from] sqlx::Error),
}

impl IntoResponse for ApiError {
    fn into_response(self) -> Response {
        let (status, code, message) = match &self {
            ApiError::BadRequest(m) => (StatusCode::BAD_REQUEST, "bad_request", m.clone()),
            ApiError::NotFound(m) => (StatusCode::NOT_FOUND, "not_found", m.clone()),
            ApiError::InvalidStructure(m) => {
                (StatusCode::UNPROCESSABLE_ENTITY, "invalid_structure", m.clone())
            }
            // Not every gRPC failure is an outage. The worker aborts with
            // INVALID_ARGUMENT when the user's SMILES cannot be parsed -- that
            // is a fact about the input, and reporting it as "service
            // unavailable" tells someone with a typo that our server is broken.
            ApiError::ChemWorker(status) => match status.code() {
                tonic::Code::InvalidArgument => (
                    StatusCode::UNPROCESSABLE_ENTITY,
                    "invalid_structure",
                    status.message().to_string(),
                ),
                tonic::Code::NotFound => (
                    StatusCode::NOT_FOUND,
                    "not_found",
                    status.message().to_string(),
                ),
                code => {
                    tracing::error!(error = %status, ?code, "chemworker call failed");
                    (
                        StatusCode::SERVICE_UNAVAILABLE,
                        "chemworker_unavailable",
                        "The cheminformatics service is not responding.".to_string(),
                    )
                }
            },
            ApiError::Database(err) => {
                tracing::error!(error = %err, "database query failed");
                (
                    StatusCode::INTERNAL_SERVER_ERROR,
                    "internal_error",
                    "An internal error occurred.".to_string(),
                )
            }
        };

        (status, Json(json!({ "error": { "code": code, "message": message } }))).into_response()
    }
}

pub type ApiResult<T> = Result<T, ApiError>;
