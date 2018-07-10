#ifndef LINEARHAM_NEWPHYLODATA_
#define LINEARHAM_NEWPHYLODATA_

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <libpll/pll.h>
#include <Eigen/Dense>
#include <pll_partition.hpp>
#include "NewData.hpp"
#include "VDJGermline.hpp"

/// @file NewPhyloData.hpp
/// @brief Header for the NewPhyloData class.

namespace linearham {


class NewPhyloData : public NewData {
 private:
  // PhyloHMM emission probability data structures
  Eigen::MatrixXi msa_;
  Eigen::MatrixXi xmsa_;
  std::vector<std::string> xmsa_labels_;
  std::vector<std::string> xmsa_seqs_;
  int xmsa_naive_ind_;
  Eigen::VectorXd xmsa_emission_;
  pll_utree_t* tree_;
  std::unique_ptr<pt::pll::Partition> partition_;

  // PhyloHMM xMSA site indices
  Eigen::VectorXi vgerm_xmsa_inds_;
  Eigen::MatrixXi vd_junction_xmsa_inds_;
  Eigen::VectorXi dgerm_xmsa_inds_;
  Eigen::MatrixXi dj_junction_xmsa_inds_;
  Eigen::VectorXi jgerm_xmsa_inds_;

  // Initialization Functions
  void InitializeMsa(const std::vector<std::string>& msa_seqs,
                     unsigned int tip_node_count, unsigned int sites,
                     const std::string& alphabet);

  void InitializeXmsaStructs(const std::string& alphabet);

  void InitializeHMMEmission(
      const std::unordered_map<std::string, GermlineGene>& ggenes) override;

  // Auxiliary Functions
  void BuildXmsa(const std::map<std::pair<int, int>, int>& xmsa_ids,
                 const std::string& alphabet);

  void FillHMMGermlineEmission(const Eigen::VectorXi& xmsa_inds_,
                               Eigen::VectorXd& emission_);

  void FillHMMJunctionEmission(const Eigen::MatrixXi& xmsa_inds_,
                               Eigen::MatrixXd& emission_);

 public:
  NewPhyloData(const std::string& flexbounds_str, const std::string& relpos_str,
               const std::unordered_map<std::string, GermlineGene>& ggenes,
               const std::string& newick_path, const std::string& fasta_path,
               const std::string& raxml_path);
  ~NewPhyloData();

  const Eigen::MatrixXi& msa() const { return msa_; };
  const Eigen::MatrixXi& xmsa() const { return xmsa_; };
  const std::vector<std::string>& xmsa_labels() const { return xmsa_labels_; };
  const std::vector<std::string>& xmsa_seqs() const { return xmsa_seqs_; };
  int xmsa_naive_ind() const { return xmsa_naive_ind_; };
  const Eigen::VectorXd& xmsa_emission() const { return xmsa_emission_; };
  const pll_utree_t& tree() const { return *tree_; };
  const pt::pll::Partition& partition() const { return *partition_; };
  const Eigen::VectorXi& vgerm_xmsa_inds() const { return vgerm_xmsa_inds_; };
  const Eigen::MatrixXi& vd_junction_xmsa_inds() const {
    return vd_junction_xmsa_inds_;
  };
  const Eigen::VectorXi& dgerm_xmsa_inds() const { return dgerm_xmsa_inds_; };
  const Eigen::MatrixXi& dj_junction_xmsa_inds() const {
    return dj_junction_xmsa_inds_;
  };
  const Eigen::VectorXi& jgerm_xmsa_inds() const { return jgerm_xmsa_inds_; };
};


typedef std::shared_ptr<NewPhyloData> NewPhyloDataPtr;


// Constructor Function

NewPhyloDataPtr ReadNewPhyloData(const std::string& csv_path,
                                 const std::string& dir_path,
                                 const std::string& newick_path,
                                 const std::string& fasta_path,
                                 const std::string& raxml_path);


// Auxiliary Functions

void StoreGermlineXmsaIndices(const std::vector<int>& naive_bases_,
                              const std::vector<int>& site_inds_,
                              std::map<std::pair<int, int>, int>& xmsa_ids,
                              Eigen::VectorXi& xmsa_inds_);

void StoreJunctionXmsaIndices(const std::vector<int>& naive_bases_,
                              const std::vector<int>& site_inds_,
                              std::pair<int, int> left_flexbounds,
                              std::pair<int, int> right_flexbounds,
                              std::map<std::pair<int, int>, int>& xmsa_ids,
                              Eigen::MatrixXi& xmsa_inds_);

void StoreXmsaIndex(std::pair<int, int> id,
                    std::map<std::pair<int, int>, int>& xmsa_ids,
                    int& xmsa_ind);


}  // namespace linearham

#endif  // LINEARHAM_NEWPHYLODATA_
