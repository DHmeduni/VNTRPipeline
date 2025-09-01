#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <filesystem>
#include <regex>
#include <limits>
#include <tuple>
#include <unordered_set>

// Convert string to uppercase
void to_uppercase(std::string& s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::toupper(c); });
}

// Levenshtein distance with early exit if distance > max_dist
int levenshtein(const std::string& s1, const std::string& s2, int max_dist) {
    size_t len1 = s1.size(), len2 = s2.size();
    if (std::abs((int)len1 - (int)len2) > max_dist) return max_dist + 1;
    std::vector<int> prev(len2 + 1), curr(len2 + 1);
    for (size_t j = 0; j <= len2; ++j) prev[j] = j;
    for (size_t i = 1; i <= len1; ++i) {
        curr[0] = i;
        int min_curr = curr[0];
        for (size_t j = 1; j <= len2; ++j) {
            int cost = (s1[i-1] == s2[j-1]) ? 0 : 1;
            curr[j] = std::min({ prev[j] + 1, curr[j-1] + 1, prev[j-1] + cost });
            min_curr = std::min(min_curr, curr[j]);
        }
        if (min_curr > max_dist) return max_dist + 1; // Early exit
        std::swap(prev, curr);
    }
    return prev[len2];
}

// Reverse complement of DNA
std::string revcomp(const std::string& seq) {
    std::string rc(seq.rbegin(), seq.rend());
    for (char& c : rc) {
        switch (std::toupper(c)) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            default: break;
        }
    }
    return rc;
}

// Decompress gzipped file to string
std::string decompress_gz_file(const std::string& filename) {
    gzFile file = gzopen(filename.c_str(), "rb");
    if (!file) throw std::runtime_error("Cannot open gzipped file: " + filename);
    std::string content;
    char buffer[4096];
    int num_read;
    while ((num_read = gzread(file, buffer, sizeof(buffer))) > 0) {
        content.append(buffer, num_read);
    }
    gzclose(file);
    return content;
}

// Parse FASTA content into vector of (header, sequence)
std::vector<std::pair<std::string, std::string>> parse_fasta(const std::string& fasta_content) {
    std::vector<std::pair<std::string, std::string>> records;
    std::istringstream iss(fasta_content);
    std::string line, header, seq;
    while (std::getline(iss, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!header.empty()) {
                records.emplace_back(header, seq);
                seq.clear();
            }
            header = line.substr(1);
        } else {
            seq += line;
        }
    }
    if (!header.empty()) {
        records.emplace_back(header, seq);
    }
    return records;
}

// Get all k-mers with offsets from motif
std::vector<std::pair<std::string, size_t>> get_kmers_with_offsets(const std::string& seq, int k) {
    std::vector<std::pair<std::string, size_t>> kmers;
    for (size_t i = 0; i + k <= seq.size(); ++i)
        kmers.emplace_back(seq.substr(i, k), i);
    return kmers;
}

// Find best Levenshtein match in window by sliding motif over all subwindows of allowed lengths, with early exit
std::tuple<int, size_t, size_t> best_levenshtein_in_window(
    const std::string& window,
    const std::string& motif,
    int max_mismatches)
{
    int min_dist = std::numeric_limits<int>::max();
    size_t best_offset = std::string::npos;
    size_t best_len = 0;
    size_t motif_len = motif.size();
    size_t min_len = (motif_len > (size_t)max_mismatches) ? motif_len - max_mismatches : 1;
    size_t max_len = motif_len + max_mismatches;

    for (size_t sub_len = min_len; sub_len <= max_len; ++sub_len) {
        if (window.size() < sub_len) continue;
        for (size_t offset = 0; offset + sub_len <= window.size(); ++offset) {
            std::string sub = window.substr(offset, sub_len);
            int dist = levenshtein(sub, motif, max_mismatches);
            if (dist < min_dist) {
                min_dist = dist;
                best_offset = offset;
                best_len = sub_len;
            }
            if (min_dist == 0) break; // Early exit for perfect match in this subwindow length
        }
        if (min_dist == 0) break; // Early exit for perfect match in any subwindow length
    }
    return {min_dist, best_offset, best_len};
}

// Struct to hold motif match info
struct MotifMatch {
    size_t pos;         // Position in sequence
    int distance;       // Levenshtein distance
    std::string match;  // Matched substring
    size_t length;      // Length of matched substring
};

// Find best fuzzy match of motif in sequence using k-mer seeding and both-side padding, with k-mer count pre-filter
MotifMatch find_best_fuzzy_kmer(
    const std::string& seq,
    const std::string& motif,
    int max_mismatches,
    int k,
    int pad,
    bool is_forward,
    int kmer_threshold = 15 // Default: at least 1 motif k-mer in window
)
{
    MotifMatch best = {std::string::npos, std::numeric_limits<int>::max(), "", 0};
    if (motif.size() < (size_t)k) k = motif.size();
    auto motif_kmers = get_kmers_with_offsets(motif, k);
    std::unordered_set<std::string> motif_kmer_set;
    for (const auto& [kmer, _] : motif_kmers) motif_kmer_set.insert(kmer);

    size_t seq_len = seq.size();
    size_t motif_len = motif.size();

    for (size_t pos = 0; pos + k <= seq_len; ++pos) {
        std::string kmer = seq.substr(pos, k);
        for (const auto& [motif_kmer, motif_kmer_offset] : motif_kmers) {
            if (kmer == motif_kmer) {
                size_t candidate_start = (pos >= motif_kmer_offset) ? pos - motif_kmer_offset : 0;
                size_t window_start = (candidate_start >= pad) ? candidate_start - pad : 0;
                size_t window_end = std::min(candidate_start + motif_len + pad, seq_len);
                size_t window_len = window_end - window_start;
                if (window_len < 1) continue;
                std::string window = seq.substr(window_start, window_len);

                // --- k-mer count pre-filter ---
                int kmer_matches = 0;
                for (size_t i = 0; i + k <= window.size(); ++i) {
                    if (motif_kmer_set.count(window.substr(i, k))) ++kmer_matches;
                }
                if (kmer_matches < kmer_threshold) continue; // Skip Levenshtein

                auto [dist, offset, best_len] = best_levenshtein_in_window(window, motif, max_mismatches);
                if (dist < best.distance) {
                    best.distance = dist;
                    best.pos = window_start + offset;
                    best.length = best_len;
                    best.match = seq.substr(best.pos, best_len);

                    std::cerr << "[DEBUG] " << (is_forward ? "FORWARD" : "REVERSE")
                              << " Window: '" << window << "' Motif: '" << motif
                              << "' Best offset: " << offset
                              << " Best length: " << best_len
                              << " Distance: " << dist
                              << " Best match: '" << best.match << "'\n";
                }
                if (best.distance == 0) return best; // Early exit for perfect match in this window
            }
        }
    }
    return best;
}

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <work_dir> [input_base] <start_seq> <end_seq> <max_mismatches>\n";
        return 1;
    }

    std::string work_dir = argv[1];
    std::vector<std::string> input_bases;
    int arg_offset = 0;

    if (argc < 3 || std::string(argv[2]).empty()) {
        if (argc < 5) {
            std::cerr << "Error: Not enough arguments. Usage: " << argv[0] << " <work_dir> [input_base] <start_seq> <end_seq> <max_mismatches>\n";
            return 1;
        }
        std::regex haplotype_regex("haplotype[12]$");
        for (const auto& entry : fs::directory_iterator(work_dir)) {
            if (entry.is_directory()) {
                std::string dirname = entry.path().filename().string();
                if (std::regex_search(dirname, haplotype_regex)) {
                    input_bases.push_back(dirname);
                }
            }
        }
        if (input_bases.empty()) {
            std::cerr << "No directory ending with 'haplotype1' or 'haplotype2' found in " << work_dir << "\n";
            return 1;
        }
        arg_offset = 1;
    } else {
        if (argc < 6) {
            std::cerr << "Error: Not enough arguments. Usage: " << argv[0] << " <work_dir> [input_base] <start_seq> <end_seq> <max_mismatches>\n";
            return 1;
        }
        input_bases.push_back(argv[2]);
    }

    std::string start_seq = argv[2 + arg_offset];
    std::string end_seq = argv[3 + arg_offset];
    int max_mismatches = std::stoi(argv[4 + arg_offset]);
    int k = 8; // default k-mer size
    int pad = max_mismatches; // both-side padding

    to_uppercase(start_seq);
    to_uppercase(end_seq);

    std::string start_seq_rev = revcomp(start_seq);
    std::string end_seq_rev = revcomp(end_seq);
    to_uppercase(start_seq_rev);
    to_uppercase(end_seq_rev);

    for (const auto& input_base : input_bases) {
        std::string input_file = work_dir + "/" + input_base + "/" + input_base + ".trimmedReads.fasta.gz";
        std::string output_for = work_dir + "/" + input_base + "/" + input_base + ".trimmedReads_for.fasta";
        std::string output_rev = work_dir + "/" + input_base + "/" + input_base + ".trimmedReads_rev.fasta";

        try {
            if (!fs::exists(input_file)) {
                throw std::runtime_error("Input file does not exist: " + input_file);
            }

            std::string file_content = decompress_gz_file(input_file);
            auto records = parse_fasta(file_content);

            std::ofstream out_for(output_for);
            std::ofstream out_rev(output_rev);

            int for_count = 0, rev_count = 0;
            for (auto& [header, seq] : records) {
                to_uppercase(seq);

                // Forward search
                MotifMatch start_match = find_best_fuzzy_kmer(seq, start_seq, max_mismatches, k, pad, true);
                if (start_match.pos != std::string::npos && start_match.distance <= max_mismatches) {
                    size_t search_from = start_match.pos + start_match.length;
                    if (search_from < seq.size()) {
                        MotifMatch end_match = find_best_fuzzy_kmer(seq.substr(search_from), end_seq, max_mismatches, k, pad, true);
                        if (end_match.pos != std::string::npos && end_match.distance <= max_mismatches) {
                            size_t end_pos = search_from + end_match.pos;
                            size_t extract_len = end_pos + end_match.length - start_match.pos;
                            if (start_match.pos + extract_len <= seq.size()) {
                                std::string region = seq.substr(start_match.pos, extract_len);
                                out_for << ">" << header << "\n" << region << "\n";
                                ++for_count;
                            }
                        }
                    }
                }

                // Reverse search
                MotifMatch end_match_rev = find_best_fuzzy_kmer(seq, end_seq_rev, max_mismatches, k, pad, false);
                if (end_match_rev.pos != std::string::npos && end_match_rev.distance <= max_mismatches) {
                    size_t search_from = end_match_rev.pos + end_match_rev.length;
                    if (search_from < seq.size()) {
                        MotifMatch start_match_rev = find_best_fuzzy_kmer(seq.substr(search_from), start_seq_rev, max_mismatches, k, pad, false);
                        if (start_match_rev.pos != std::string::npos && start_match_rev.distance <= max_mismatches) {
                            size_t start_pos = search_from + start_match_rev.pos;
                            size_t extract_len = start_pos + start_match_rev.length - end_match_rev.pos;
                            if (end_match_rev.pos + extract_len <= seq.size()) {
                                std::string region = seq.substr(end_match_rev.pos, extract_len);
                                out_rev << ">" << header << "\n" << region << "\n";
                                ++rev_count;
                            }
                        }
                    }
                }
            }
            out_for.close();
            out_rev.close();

            std::cout << "Done for " << input_base << ". Forward: " << for_count << " Reverse: " << rev_count << std::endl;
        } catch (const std::exception& ex) {
            std::cerr << "Error processing " << input_base << ": " << ex.what() << std::endl;
        }
    }
    return 0;
}
