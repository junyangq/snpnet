// [[Rcpp::depends(BH)]]

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/exceptions.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <Rcpp.h>
#include <string>

class BEDMatrixPlus {
    public:
        BEDMatrixPlus(std::string path, std::size_t n, std::size_t p);
        Rcpp::NumericMatrix multiply_residuals(Rcpp::IntegerVector i, Rcpp::IntegerVector j, Rcpp::NumericVector imputeValues, Rcpp::NumericMatrix residuals);
    private:
        BEDMatrixPlus(const BEDMatrixPlus&);
        BEDMatrixPlus& operator=(const BEDMatrixPlus&);
        boost::interprocess::file_mapping file;
        boost::interprocess::mapped_region file_region;
        const char* file_data;
        std::size_t nrow;
        std::size_t ncol;
        unsigned short int byte_padding; // Each new "row" starts a new byte
        static const unsigned short int length_header;
};


BEDMatrixPlus::BEDMatrixPlus(std::string path, std::size_t n, std::size_t p) : nrow(n), ncol(p), byte_padding((n % 4 == 0) ? 0 : 4 - (n % 4)) {
    try {
        this->file = boost::interprocess::file_mapping(path.c_str(), boost::interprocess::read_only);
    } catch(const boost::interprocess::interprocess_exception& e) {
        throw std::runtime_error("File not found.");
    }
    this->file_region = boost::interprocess::mapped_region(this->file, boost::interprocess::read_only);
    this->file_data = static_cast<const char*>(this->file_region.get_address());
    // Check magic number
    if (!(this->file_data[0] == '\x6C' && this->file_data[1] == '\x1B')) {
        throw std::runtime_error("File is not a binary PED file.");
    }
    // Check mode: 00000001 indicates the default variant-major mode (i.e.
    // list all samples for first variant, all samples for second variant,
    // etc), 00000000 indicates the unsupported sample-major mode (i.e. list
    // all variants for the first sample, list all variants for the second
    // sample, etc)
    if (this->file_data[2] != '\x01') {
        throw std::runtime_error("Sample-major mode is not supported.");
    }
    // Get number of bytes
    const std::size_t num_bytes = this->file_region.get_size();
    // Check if given dimensions match the file

    // ADDDDDDD
    if ((this->nrow * this->ncol) + (this->byte_padding * this->ncol) != (num_bytes - this->length_header) * 4) {
        throw std::runtime_error("n or p does not match the dimensions of the file.");
    }
}


Rcpp::NumericMatrix BEDMatrixPlus::multiply_residuals(Rcpp::IntegerVector i, Rcpp::IntegerVector j, Rcpp::NumericVector imputeValues, Rcpp::NumericMatrix residuals) {

    if (Rcpp::is_true(Rcpp::any(i > this->nrow)) || Rcpp::is_true(Rcpp::any(j > this->ncol))) {
        throw std::runtime_error("subscript out of bounds");
    }

    Rcpp::IntegerVector i0(i - 1);
    Rcpp::IntegerVector j0(j - 1);

    std::size_t size_i = i.size();
    std::size_t size_j = j.size();
    std::size_t size_k = residuals.ncol();

    Rcpp::NumericMatrix out(size_j, size_k);

    for (std::size_t idx_j = 0; idx_j < size_j; idx_j++) {

        double summation[size_k][4];
        // double doubleType;
        std::memset(*summation, 0, size_k * 4 * sizeof(double));

        std::size_t base_pos = (j0[idx_j] * this->nrow) + (this->byte_padding * j0[idx_j]);

        for (std::size_t idx_i = 0; idx_i < size_i; idx_i++) {

            std::size_t which_pos = base_pos + i0[idx_i];

            std::size_t which_byte = std::floor(which_pos / 4);
                // Find genotype in byte
            unsigned short int which_genotype = (which_pos % 4) * 2;
                // Read in the whole byte
            char genotypes = this->file_data[which_byte + this->length_header];
                // Remove the other genotypes by shifting the genotype of interest
                // to the end of the byte and masking with 00000011
            char genotype = genotypes >> which_genotype & 3;
                // Remap genotype value to resemble RAW file, i.e. 0 indicates homozygous
                // major allele, 1 indicates heterozygous, and 2 indicates homozygous minor
                // allele. In BED, the coding is different: homozygous minor allele is 0
                // (00) and homozygous major allele is 3 (11). Each byte is read backwards.

            if (genotype != 3) {
                for (std::size_t idx_k = 0; idx_k < size_k; idx_k++) {
                    summation[idx_k][genotype] += residuals(i0[idx_i], idx_k);
                }
            }
        }

        for (std::size_t idx_k = 0; idx_k < size_k; idx_k++) {
            out(idx_j, idx_k) = summation[idx_k][2] + 2 * summation[idx_k][0] + imputeValues(j0[idx_j]) * summation[idx_k][1];
        }
    }
    return out;
}

const unsigned short int BEDMatrixPlus::length_header = 3;

RcppExport SEXP BEDMatrix__multiply_residuals(SEXP xp_, SEXP n_, SEXP p_, SEXP i_, SEXP j_, SEXP iv_, SEXP r_) {
    std::string path = Rcpp::as<std::string>(xp_);
    std::size_t n = Rcpp::as<std::size_t>(n_);
    std::size_t p = Rcpp::as<std::size_t>(p_);
    Rcpp::IntegerVector i = Rcpp::as<Rcpp::IntegerVector>(i_);
    Rcpp::IntegerVector j = Rcpp::as<Rcpp::IntegerVector>(j_);
    Rcpp::NumericVector iv = Rcpp::as<Rcpp::NumericVector>(iv_);
    Rcpp::NumericMatrix r = Rcpp::as<Rcpp::NumericMatrix>(r_);
    try {
        Rcpp::XPtr<BEDMatrixPlus> ptr(new BEDMatrixPlus(path, n, p), true);
        Rcpp::NumericMatrix res = ptr->multiply_residuals(i, j, iv, r);
        return res;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
        return 0;
    }
};