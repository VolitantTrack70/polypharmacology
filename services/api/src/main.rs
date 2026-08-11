//! Polypharmacology knowledge graph -- HTTP API.
//!
//! This service orchestrates; it does not compute. Structure handling is
//! delegated to the Python chemworker over gRPC (RDKit has no mature Rust
//! equivalent), and the relational/graph queries go to Postgres. What lives
//! here is request validation, the two-phase query, and response assembly.

mod chemworker;
mod error;
mod routes;

use std::time::Duration;

use axum::routing::{get, post};
use axum::Router;
use sqlx::postgres::PgPoolOptions;
use sqlx::PgPool;
use tower_http::compression::CompressionLayer;
use tower_http::cors::CorsLayer;
use tower_http::trace::TraceLayer;
use tracing_subscriber::{layer::SubscriberExt, util::SubscriberInitExt, EnvFilter};

use crate::chemworker::ChemWorker;

#[derive(Clone)]
pub struct AppState {
    pub db: PgPool,
    pub chem: ChemWorker,
    pub defaults: Defaults,
}

/// Defaults applied when a request omits them.
#[derive(Clone, Copy)]
pub struct Defaults {
    /// See docs/decisions/0002 -- 0.40, deliberately not 0.85.
    pub tanimoto: f32,
    pub pchembl: f32,
    pub similar_limit: u32,
}

impl Default for Defaults {
    fn default() -> Self {
        Self { tanimoto: 0.40, pchembl: 6.0, similar_limit: 250 }
    }
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    dotenvy::from_filename("../../.env").ok();

    tracing_subscriber::registry()
        .with(EnvFilter::try_from_default_env()
            .unwrap_or_else(|_| "info,chemmed_api=debug".into()))
        .with(tracing_subscriber::fmt::layer())
        .init();

    let database_url = std::env::var("DATABASE_URL")
        .unwrap_or_else(|_| "postgres://chemmed:change_me_locally@localhost:5432/chemmed".into());
    let chem_addr = std::env::var("CHEMWORKER_GRPC_ADDR")
        .unwrap_or_else(|_| "http://localhost:50051".into());
    let bind_addr = std::env::var("API_BIND_ADDR")
        .unwrap_or_else(|_| "0.0.0.0:8080".into());

    let db = PgPoolOptions::new()
        .max_connections(16)
        .acquire_timeout(Duration::from_secs(10))
        .connect(&database_url)
        .await?;
    tracing::info!("connected to postgres");

    // Lazy: the worker loads a ~600MB index and may still be starting.
    let chem = ChemWorker::connect_lazy(&chem_addr)?;
    tracing::info!(addr = %chem_addr, "chemworker client ready (lazy)");

    let state = AppState { db, chem, defaults: Defaults::default() };

    let api = Router::new()
        .route("/health", get(routes::health::health))
        .route("/status", get(routes::health::status))
        .route("/resolve", post(routes::compound::resolve))
        .route("/compound/{chembl_id}", get(routes::compound::get_compound))
        .route("/offtargets", post(routes::offtargets::off_targets));

    let app = Router::new()
        .nest("/api", api)
        .layer(TraceLayer::new_for_http())
        .layer(CompressionLayer::new())
        // Dev-permissive. Tighten before this is reachable off-machine.
        .layer(CorsLayer::permissive())
        .with_state(state);

    let listener = tokio::net::TcpListener::bind(&bind_addr).await?;
    tracing::info!("listening on http://{bind_addr}");

    axum::serve(listener, app)
        .with_graceful_shutdown(shutdown_signal())
        .await?;

    Ok(())
}

async fn shutdown_signal() {
    let _ = tokio::signal::ctrl_c().await;
    tracing::info!("shutting down");
}
