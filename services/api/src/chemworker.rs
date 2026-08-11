//! gRPC client for the Python cheminformatics worker.

use tonic::transport::Channel;

pub mod pb {
    //! Generated from `services/chemworker/proto/chemworker.proto`.
    tonic::include_proto!("chemmed.chemworker.v1");
}

use pb::chem_worker_client::ChemWorkerClient;
use pb::{
    similarity_search_request::Query, SimilaritySearchRequest, SimilaritySearchResponse,
    StandardizeRequest, StandardizeResponse, StatusRequest, StatusResponse,
};

/// Thin wrapper over the generated client.
///
/// `ChemWorkerClient<Channel>` is cheap to clone -- the underlying channel is
/// reference-counted and multiplexes over one HTTP/2 connection -- so handlers
/// clone this out of shared state rather than reconnecting per request.
#[derive(Clone)]
pub struct ChemWorker {
    client: ChemWorkerClient<Channel>,
}

impl ChemWorker {
    /// Connects lazily. The worker holds a ~600MB fingerprint index and can
    /// take a while to come up, so the API must not require it at boot.
    pub fn connect_lazy(addr: &str) -> anyhow::Result<Self> {
        let endpoint = Channel::from_shared(addr.to_string())
            .map_err(|e| anyhow::anyhow!("invalid chemworker address {addr:?}: {e}"))?;
        Ok(Self {
            client: ChemWorkerClient::new(endpoint.connect_lazy()),
        })
    }

    pub async fn standardize(&self, smiles: &str) -> Result<StandardizeResponse, tonic::Status> {
        let mut client = self.client.clone();
        let resp = client
            .standardize(StandardizeRequest {
                smiles: smiles.to_string(),
            })
            .await?;
        Ok(resp.into_inner())
    }

    pub async fn search_by_smiles(
        &self,
        smiles: &str,
        threshold: f32,
        limit: u32,
    ) -> Result<SimilaritySearchResponse, tonic::Status> {
        self.search(Query::Smiles(smiles.to_string()), threshold, limit)
            .await
    }

    pub async fn search_by_chembl_id(
        &self,
        chembl_id: &str,
        threshold: f32,
        limit: u32,
    ) -> Result<SimilaritySearchResponse, tonic::Status> {
        self.search(Query::ChemblId(chembl_id.to_string()), threshold, limit)
            .await
    }

    async fn search(
        &self,
        query: Query,
        threshold: f32,
        limit: u32,
    ) -> Result<SimilaritySearchResponse, tonic::Status> {
        let mut client = self.client.clone();
        let resp = client
            .similarity_search(SimilaritySearchRequest {
                query: Some(query),
                threshold,
                limit,
            })
            .await?;
        Ok(resp.into_inner())
    }

    pub async fn status(&self) -> Result<StatusResponse, tonic::Status> {
        let mut client = self.client.clone();
        Ok(client.status(StatusRequest {}).await?.into_inner())
    }
}
