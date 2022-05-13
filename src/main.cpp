// Copyright (c) 2022 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <string_view>
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <tuple>
#include <iostream>
#include <memory>
#include <optional>
#include <filesystem>
#include <variant>
#include <type_traits>
#include <array>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "version.hpp"

namespace fs = std::filesystem;

struct HtsFileDeleter
{
    void operator()(htsFile* file) const { hts_close(file); }
};
struct HtsHeaderDeleter
{
    void operator()(bam_hdr_t* header) const { bam_hdr_destroy(header); }
};
struct HtsIndexDeleter
{
    void operator()(hts_idx_t* index) const { hts_idx_destroy(index); }
};
struct HtsIteratorDeleter
{
    void operator()(hts_itr_t* iterator) const { sam_itr_destroy(iterator); }
};
struct HtsBam1Deleter
{
    void operator()(bam1_t* b) const { bam_destroy1(b); }
};

using ReadNameMap = std::unordered_map<std::string, std::string>;

auto estimate_line_count(const fs::path& filename)
{
    const auto bytes = fs::file_size(filename);
    std::ifstream file {filename};
    std::string line {};
    std::getline(file, line);
    return bytes / (line.length() + 1);
}

auto load_reads(const fs::path& qnames_tsv_path, const bool verbose = false)
{
    const std::size_t log_tick {10'000'000};
    const auto estimated_lines = estimate_line_count(qnames_tsv_path);
    std::ifstream qname_tsv {qnames_tsv_path};
    std::vector<std::pair<std::string, std::string>> reads {};
    reads.reserve(estimated_lines);
    std::string line {};
    for (std::size_t i {0}; std::getline(qname_tsv, line); ++i) {
        if (verbose && i > 0 and i % log_tick == 0) {
            std::clog << "Loaded " << i << " reads" << std::endl;
        }
        const auto name_end_pos = line.find('\t');
        reads.emplace_back(line.substr(0, name_end_pos), name_end_pos != std::string::npos ? line.substr(name_end_pos + 1) : "");
    }
    return ReadNameMap {
        std::make_move_iterator(std::begin(reads)),
        std::make_move_iterator(std::end(reads))
    };
}

struct Tag
{
    using Name = std::array<char, 2>;
    using Value = std::variant<std::string_view, long, float>;
    Name id;
    Value value;
};

Tag::Value get_tag_value(const std::string_view value)
{
    try {
        if (value.find('.') == std::string::npos) {
            return std::stol(std::string{value});
        } else {
            return std::stof(std::string{value});
        }
    } catch (const std::invalid_argument&) {
        return value;
    }
}

auto get_tag(const std::string_view& tag)
{
    if (tag.size() < 2 || tag.size() == 3 || (tag.size() > 3 && tag[2] != ':')) {
        std::clog << "Invalid tag " << tag << " (required TAG:VALUE)" << std::endl;
        exit(1);
    }
    Tag result {};
    std::copy_n(std::cbegin(tag), 2, std::begin(result.id));
    if (tag.size() > 3) {
        result.value = get_tag_value(tag.substr(3));
    }
    return result;
}

template<class> inline constexpr bool always_false_v = false;

void add_tag(const Tag& tag, bam1_t* rec)
{
    std::visit([&] (auto&& value) {
        using T = std::remove_cv_t<std::decay_t<decltype(value)>>;
        if constexpr (std::is_same_v<T, long>) {
            bam_aux_update_int(rec, tag.id.data(), value);
        } else if constexpr (std::is_same_v<T, float>) {
            bam_aux_update_float(rec, tag.id.data(), value);
        } else if constexpr (std::is_same_v<T, std::string_view>) {
            if (!value.empty()) {
                bam_aux_update_str(rec, tag.id.data(), value.size() + 1, value.data());
            }
        } else { 
            static_assert(always_false_v<T>, "non-exhaustive visitor!");
        }
    }, tag.value);
}

void 
tag_reads(const fs::path& src_bam_path,
          const ReadNameMap& read_names,
          const std::optional<fs::path>& dst_bam_path,
          const std::optional<Tag>& tag = std::nullopt,
          const std::optional<std::uint16_t> flag = std::nullopt,
          const bool verbose = false)
{
    const std::size_t log_tick {10'000'000};
    std::unique_ptr<htsFile, HtsFileDeleter> src_bam {sam_open(src_bam_path.c_str(), "r"), HtsFileDeleter {}};
    std::unique_ptr<bam_hdr_t, HtsHeaderDeleter> header {sam_hdr_read(src_bam.get()), HtsHeaderDeleter {}};
    std::unique_ptr<htsFile, HtsFileDeleter> dst_bam {
        dst_bam_path ? sam_open(dst_bam_path->c_str(), "[w]b") : sam_open("-", "w"), 
        HtsFileDeleter {}};
    if (sam_hdr_write(dst_bam.get(), header.get()) < 0) {
        std::clog << "Error writing BAM" << std::endl;
        exit(1);
    }
    std::unique_ptr<bam1_t, HtsBam1Deleter> rec {bam_init1(), HtsBam1Deleter {}};
    int r;
    std::string name {};
    std::size_t i {0}, marked {0};
    std::vector<Tag> tags {};
    bool set_default_tag_value {};
    if (tag) {
        tags.push_back(*tag);
        try {
            // If no value is provided for the default tag then
            // the input file contains values
            set_default_tag_value = std::get<std::string_view>(tag->value).empty();
        } catch (const std::bad_variant_access&) {}
    }
    while ((r = sam_read1(src_bam.get(), header.get(), rec.get())) >= 0) {
        if (verbose && i > 0 and i % log_tick == 0) {
            std::clog << "Processed " << i << " reads -- marked " << marked
                      << " (~" << static_cast<int>(100 * marked / i) << "%)"
                      << std::endl;
        }
        name = bam_get_qname(rec.get());
        const auto read_itr = read_names.find(name);
        if (read_itr != std::cend(read_names)) {
            const std::string& edits {read_itr->second};
            auto read_flag = flag;
            if (!edits.empty()) {
                const auto tags_end = edits.find('\t');
                std::string_view reads_tags {edits};
                if (tags_end != std::string::npos) {
                    auto f = static_cast<std::uint16_t>(std::stoi(edits.substr(tags_end + 1)));
                    if (read_flag) {
                        *read_flag |= f;
                    } else {
                        read_flag = f;
                    }
                    reads_tags = edits.substr(0, tags_end);
                }
                // TODO - accept multiple tags seperated with ;
                if (set_default_tag_value) {
                    tags.front().value = get_tag_value(reads_tags);
                } else {
                    tags.push_back(get_tag(reads_tags));
                }
            }
            if (read_flag) {
                rec->core.flag |= *read_flag;
            }
            if (tags.empty() && !read_flag) {
                std::clog << "WARN: no tags or flags for read " << name << std::endl;
            }
            for (const auto& t : tags) {
                add_tag(t, rec.get());
            }
            if (tag) {
                tags.resize(1);
            } else {
                tags.clear();
            }
            ++marked;
        }
        if (sam_write1(dst_bam.get(), header.get(), rec.get()) < 0) {
            std::clog << "Error writing BAM" << std::endl;
            exit(1);
        }
        ++i;
    }
    if (verbose) {
        std::clog << "Processed " << i << " reads -- marked " << marked
                  << " (~" << static_cast<int>(100 * marked / i) << "%)"
                  << std::endl;
    }
}

struct ProgramOptions
{
    fs::path src_bam_path, qname_tsv_path;
    std::optional<fs::path> output;
    std::optional<Tag> tag;
    std::optional<std::uint16_t> flag;
    bool build_index;
    int verbose = 0;
};

void print_usage(const std::string_view& program)
{
    std::cout << "Usage: " << program << " [options] <in.bam> <qnames.tsv>" << std::endl;
}
void print_help()
{
    std::cout << "Options:" << std::endl;
    std::cout << "  --help              Print help information" << std::endl;
    std::cout << "  -o, --output FILE   Output bam/cram FILE" << std::endl;
    std::cout << "  -t, --tag STR1:STR2 Add tag STR1 with value STR2 to selected reads" << std::endl;
    std::cout << "  -f, --flag FLAG     Add FLAG to selected reads" << std::endl;
    std::cout << "  -i, --index         Build the bam/cram index for the output" << std::endl;
    std::cout << "  --verbose INT       Print verbose logging info [0]" << std::endl;
    std::cout << "  --version           Print version information";
}
void print_version(const std::string_view& program)
{
    std::cout << program << " " << VERSION_MAJOR << "." << VERSION_MINOR;
}

auto parse_args(const int argc, char** argv)
{
    std::string_view program {argv[0]};
    if (program.starts_with("./")) {
        program.remove_prefix(2);
    }
    constexpr int num_positionals {2};
    ProgramOptions result {};
    std::vector<std::string_view> args(argv + 1, argv + argc), positionals {};
    std::optional<std::string_view> option {};
    for (int consumed {1}; const auto& arg : args) {
        if (arg == "-h" || arg == "--help") {
            print_usage(program);
            std::cout << std::endl;
            print_help();
            std::cout << std::endl;
            exit(0);
        } else if (arg == "--version") {
            print_version(program);
            std::cout << std::endl;
            exit(0);
        } else if (argc - consumed <= num_positionals - static_cast<int>(positionals.size())) {
            positionals.push_back(arg);
        } else if (option == "-o" || option == "--output") {
            result.output = arg;
            option = std::nullopt;
        } else if (option == "-t" || option == "--tag") {
            result.tag = get_tag(arg);
            option = std::nullopt;
        } else if (option == "-f" || option == "--flag") {
            result.flag = static_cast<std::uint16_t>(std::stoi(std::string(arg)));
            option = std::nullopt;
        } else if (option == "--verbose") {
            result.verbose = std::stoi(std::string(arg));
            option = std::nullopt;
        } else if (arg == "-i" || arg == "--index") {
            result.build_index = true;
        } else if (arg.starts_with("-")) {
            option = arg;
        } else {
            positionals.push_back(arg);
        }
        ++consumed;
    }
    if (positionals.size() != num_positionals) {
        print_usage(program);
        exit(1);
    }
    result.src_bam_path = positionals[0];
    result.qname_tsv_path = positionals[1];
    return result;
}

void check_input(const ProgramOptions& options)
{
    if (options.src_bam_path != "-" && !fs::exists(options.src_bam_path)) {
        std::cerr << "ERROR: input file " << options.src_bam_path << " does not exist." << std::endl;
        exit(1);
    }
    if (!fs::exists(options.qname_tsv_path)) {
        std::cerr << "ERROR: input file " << options.qname_tsv_path << " does not exist." << std::endl;
        exit(1);
    }
    if (options.build_index && !options.output) {
        std::clog << "Warn: cannot build bam index without --output!" << std::endl;
    }
}

int main(int argc, char** argv)
{
    const auto options = parse_args(argc, argv);
    check_input(options);
    auto read_names = load_reads(options.qname_tsv_path, options.verbose);
    if (options.verbose > 0) std::clog << "Loaded " << read_names.size() << " read names" << std::endl;
    tag_reads(options.src_bam_path, read_names, options.output, options.tag, options.flag, options.verbose);
    if (options.build_index && options.output && sam_index_build(options.output->c_str(), 0) < 0) {
        std::clog << "Failed to build bam index" << std::endl;
    }
}
