// [[Rcpp::depends(BH)]]

#define PLINK_BED_HEADER_LENGTH 3
#define PLINK_BED_GENOTYPES_PER_BYTE 4

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
        Rcpp::NumericMatrix extract_matrix(Rcpp::IntegerVector i, Rcpp::IntegerVector j);
        Rcpp::NumericMatrix extract_columns(Rcpp::IntegerVector j);
        Rcpp::NumericMatrix multiply_residuals(Rcpp::IntegerVector js, Rcpp::IntegerVector je, Rcpp::NumericVector imputeValues, Rcpp::NumericMatrix residuals);
    private:
        BEDMatrixPlus(const BEDMatrixPlus&);
        BEDMatrixPlus& operator=(const BEDMatrixPlus&);
        boost::interprocess::file_mapping file;
        boost::interprocess::mapped_region file_region;
        uint8_t* file_data;
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
    this->file_data = static_cast<uint8_t*>(this->file_region.get_address());
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


Rcpp::NumericMatrix BEDMatrixPlus::multiply_residuals(Rcpp::IntegerVector js, Rcpp::IntegerVector je, Rcpp::NumericVector imputeValues, Rcpp::NumericMatrix residuals) {

    if (Rcpp::is_true(Rcpp::any(js > this->ncol)) || Rcpp::is_true(Rcpp::any(je > this->ncol))) {
        throw std::runtime_error("subscript out of bounds");
    }

    // Rcpp::IntegerVector i0(i - 1);
    // Rcpp::IntegerVector j0(j - 1);
    
    std::size_t js0 = js[0] - 1;
    std::size_t je0 = je[0] - 1;
    std::size_t size_j = je0 - js0 + 1;
    std::size_t size_i = residuals.nrow();
    std::size_t size_k = residuals.ncol();

    Rcpp::NumericMatrix out(size_j, size_k);

    for (std::size_t idx_j = js0; idx_j <= je0; idx_j++) {

        double summation[size_k][4];
        // double doubleType;
        std::memset(*summation, 0, size_k * 4 * sizeof(double));

        std::size_t base_pos = (idx_j * this->nrow) + (this->byte_padding * idx_j);

        std::size_t which_byte = base_pos / 4;
        unsigned short int which_genotype = base_pos % 4;
        uint8_t *pgeno = this->file_data + which_byte + this->length_header;
        uint8_t genotypes = *pgeno;
        genotypes >>= (which_genotype * 2);

        for (std::size_t idx_i = 0; idx_i < size_i; idx_i++) {

            // std::size_t which_byte = which_pos / 4;
                // Find genotype in byte
            // unsigned short int which_genotype = (which_pos % 4) * 2;
                // Read in the whole byte
            // char genotypes = this->file_data[which_byte + this->length_header];
                // Remove the other genotypes by shifting the genotype of interest
                // to the end of the byte and masking with 00000011
            uint8_t genotype = genotypes & 3;
                // Remap genotype value to resemble RAW file, i.e. 0 indicates homozygous
                // major allele, 1 indicates heterozygous, and 2 indicates homozygous minor
                // allele. In BED, the coding is different: homozygous minor allele is 0
                // (00) and homozygous major allele is 3 (11). Each byte is read backwards.
            if (genotype != 3) {
                for (std::size_t idx_k = 0; idx_k < size_k; idx_k++) {
                    summation[idx_k][genotype] += residuals(idx_i, idx_k);
                }
            }
            if (which_genotype != 3) {
                genotypes >>= 2;
                which_genotype += 1;
            } else {
                pgeno += 1;
                genotypes = *pgeno;
                which_genotype = 0;
            }
        }

        for (std::size_t idx_k = 0; idx_k < size_k; idx_k++) {
            out(idx_j-js0, idx_k) = summation[idx_k][2] + 2 * summation[idx_k][0] + imputeValues(idx_j) * summation[idx_k][1];
        }
    }
    return out;
}

const unsigned short int BEDMatrixPlus::length_header = 3;


RcppExport SEXP BEDMatrixPlus__new(SEXP path_, SEXP n_, SEXP p_) {
    // Convert inputs to appropriate C++ types
    std::string path = Rcpp::as<std::string>(path_);
    std::size_t n = Rcpp::as<std::size_t>(n_);
    std::size_t p = Rcpp::as<std::size_t>(p_);
    try {
        // Create a pointer to a BEDMatrix object and wrap it as an external
        // pointer
        Rcpp::XPtr<BEDMatrixPlus> ptr(new BEDMatrixPlus(path, n, p), true);
        // Return the external pointer to the R side
        return ptr;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
        return 0;
    }
}

Rcpp::NumericMatrix BEDMatrixPlus::extract_matrix(Rcpp::IntegerVector i, Rcpp::IntegerVector j) {
    // Check if indexes are out of bounds
    if (Rcpp::is_true(Rcpp::any(i > this->nrow)) && Rcpp::is_true(Rcpp::any(j > this->ncol))) {
        throw std::runtime_error("subscript out of bounds");
    }
    // Convert from 1-index to 0-index
    Rcpp::IntegerVector i0(i - 1);
    Rcpp::IntegerVector j0(j - 1);
    // Keep sizes of and j
    std::size_t size_i = i.size();
    std::size_t size_j = j.size();
    // Reserve output matrix
    Rcpp::NumericMatrix out(size_i, size_j);
    std::size_t n_bytes = this->nrow / PLINK_BED_GENOTYPES_PER_BYTE + (this->nrow % PLINK_BED_GENOTYPES_PER_BYTE != 0); // fast ceil for int

    // Iterate over column indexes
    for (std::size_t idx_j = 0; idx_j < size_j; idx_j++) {
        // Iterate over row indexes
        std::size_t prev_bytes = j0[idx_j] * n_bytes;
        for (std::size_t idx_i = 0; idx_i < size_i; idx_i++) {
            std::size_t i_bytes = i0[idx_i] / PLINK_BED_GENOTYPES_PER_BYTE;
            std::size_t i_genotypes = 2 * (i0[idx_i] - i_bytes * PLINK_BED_GENOTYPES_PER_BYTE);
            uint8_t genotypes = this->file_data[PLINK_BED_HEADER_LENGTH + (prev_bytes + i_bytes)];
            uint8_t genotype = genotypes >> i_genotypes & 3;
            double mapping = NA_REAL;
            // Load byte from map
            // Extract genotypes from byte by shifting the genotype of interest to the
            // end of the byte and masking with 00000011
            // Remap genotype value to resemble RAW file, i.e. 0 indicates homozygous
            // major allele, 1 indicates heterozygous, and 2 indicates homozygous minor
            // allele. In BED, the coding is different: homozygous minor allele is 0
            // (00) and homozygous major allele is 3 (11). Each byte is read backwards.
            // missing
            if (genotype == 0) {
                mapping = 2; // homozygous AA
            } else if (genotype == 3) {
                mapping = 0; // homozygous BB
            } else if (genotype == 2) {
                mapping = 1; // heterozygous AB
            }
            out(idx_i, idx_j) = mapping;
        }
    }
    return out;
}

Rcpp::NumericMatrix BEDMatrixPlus::extract_columns(Rcpp::IntegerVector j) {
    // Check if indexes are out of bounds
    if (Rcpp::is_true(Rcpp::any(j > this->ncol))) {
        throw std::runtime_error("subscript out of bounds");
    }
    // Convert from 1-index to 0-index
    Rcpp::IntegerVector j0(j - 1);
    // Keep sizes of and j
    std::size_t size_i = this->nrow;
    std::size_t size_j = j.size();
    // Reserve output matrix
    Rcpp::NumericMatrix out(size_i, size_j);
    std::size_t n_bytes = this->nrow / PLINK_BED_GENOTYPES_PER_BYTE + (this->nrow % PLINK_BED_GENOTYPES_PER_BYTE != 0); // fast ceil for int

    // Iterate over column indexes
    for (std::size_t idx_j = 0; idx_j < size_j; idx_j++) {
        // Iterate over row indexes
        std::size_t prev_bytes = j0[idx_j] * n_bytes;
        for (std::size_t idx_i = 0; idx_i < size_i; idx_i++) {
            std::size_t i_bytes = idx_i / PLINK_BED_GENOTYPES_PER_BYTE;
            std::size_t i_genotypes = 2 * (idx_i - i_bytes * PLINK_BED_GENOTYPES_PER_BYTE);
            uint8_t genotypes = this->file_data[PLINK_BED_HEADER_LENGTH + (prev_bytes + i_bytes)];
            uint8_t genotype = genotypes >> i_genotypes & 3;
            double mapping = NA_REAL;
            // Load byte from map
            // Extract genotypes from byte by shifting the genotype of interest to the
            // end of the byte and masking with 00000011
            // Remap genotype value to resemble RAW file, i.e. 0 indicates homozygous
            // major allele, 1 indicates heterozygous, and 2 indicates homozygous minor
            // allele. In BED, the coding is different: homozygous minor allele is 0
            // (00) and homozygous major allele is 3 (11). Each byte is read backwards.
            // missing
            if (genotype == 0) {
                mapping = 2; // homozygous AA
            } else if (genotype == 3) {
                mapping = 0; // homozygous BB
            } else if (genotype == 2) {
                mapping = 1; // heterozygous AB
            }
            out(idx_i, idx_j) = mapping;
        }
    }
    return out;
}


RcppExport SEXP BEDMatrixPlus__multiply_residuals(SEXP xp_, SEXP js_, SEXP je_, SEXP iv_, SEXP r_) {
    Rcpp::XPtr<BEDMatrixPlus> ptr(xp_);
    Rcpp::IntegerVector js = Rcpp::as<Rcpp::IntegerVector>(js_);
    Rcpp::IntegerVector je = Rcpp::as<Rcpp::IntegerVector>(je_);
    Rcpp::NumericVector iv = Rcpp::as<Rcpp::NumericVector>(iv_);
    Rcpp::NumericMatrix r = Rcpp::as<Rcpp::NumericMatrix>(r_);
    try {
        Rcpp::NumericMatrix res = ptr->multiply_residuals(js, je, iv, r);
        return res;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
        return 0;
    }
}


RcppExport SEXP BEDMatrixPlus__extract_matrix(SEXP xp_, SEXP i_, SEXP j_) {
    // Convert inputs to appropriate C++ types
    Rcpp::XPtr<BEDMatrixPlus> ptr(xp_);
    Rcpp::IntegerVector i = Rcpp::as<Rcpp::IntegerVector>(i_);
    Rcpp::IntegerVector j = Rcpp::as<Rcpp::IntegerVector>(j_);
    try {
        // Invoke the extract_matrix function
        Rcpp::NumericMatrix res = ptr->extract_matrix(i, j);
        return res;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
        return 0;
    }
}


RcppExport SEXP BEDMatrixPlus__extract_columns(SEXP xp_, SEXP j_) {
    // Convert inputs to appropriate C++ types
    Rcpp::XPtr<BEDMatrixPlus> ptr(xp_);
    Rcpp::IntegerVector j = Rcpp::as<Rcpp::IntegerVector>(j_);
    try {
        // Invoke the extract_matrix function
        Rcpp::NumericMatrix res = ptr->extract_columns(j);
        return res;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
        return 0;
    }
}
