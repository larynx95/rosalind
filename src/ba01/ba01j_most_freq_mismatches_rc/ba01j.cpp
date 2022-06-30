/*
Rosalind: BA1J
Find Frequent Words with Mismatches and Reverse Complements

We now extend “Find the Most Frequent Words with Mismatches in a String” to find
frequent words with both mismatches and reverse complements. Recall that Pattern
refers to the reverse complement of Pattern.

Frequent Words with Mismatches and Reverse Complements Problem
Find the most frequent k-mers (with mismatches and reverse complements) in a DNA
string.

Given: A DNA string Text as well as integers k and d.

Return: All k-mers Pattern maximizing the sum Countd(Text, Pattern) +
Countd(Text, Pattern) over all possible k-mers.

Sample Dataset
ACGTTGCATGTCGCATGATGCATGAGAGCT
4 1

Sample Output
ATGT ACAT
*/

#include "/home/wsl/rosalind/include/utils.hpp"

/************************************************
prototypes
************************************************/

std::string get_complement(std::string&);
int get_hamming_distance(std::string str1, std::string str2);
std::vector<std::string> perm_with_repetition(std::vector<std::string> vec,
    std::string nuc, int num);
std::vector<int> find_all_approx_occur(std::string pattern, std::string genome,
    int num_d);
std::vector<std::string> most_freq_words_with_mismatch(std::string genome,
    int num_k, int num_d);
std::vector<std::string> most_freq_words_with_mismatch_rc(std::string genome,
    int num_k, int num_d);

/************************************************
main
************************************************/

int main(int argc, const char* argv[]) {
    const std::string filename = argc < 2 ? "/home/wsl/rosalind/data/ba01j.txt" : argv[1];
    std::vector<std::string> lines = read_lines(filename);
    std::string pattern = lines[0];
    std::vector<std::string> str_nums = split_str(lines[1], " ");
    int num_k = std::stoi(str_nums[0]);
    int num_d = std::stoi(str_nums[1]);

    std::vector<std::string> answer =
        most_freq_words_with_mismatch_rc(pattern, num_k, num_d);
    print_vec(answer);

    return 0;
}

/************************************************
definitions
************************************************/

void rev_str_by_ref(std::string& str) {
    char temp;
    int len = str.length();
    for (int i = 0; i < len / 2; i++) {
        temp = str[i];
        str[i] = str[len - i - 1];
        str[len - i - 1] = temp;
    }
}

std::string rev_str(std::string& str) {
    std::string new_str;
    for (auto riter = str.rbegin(); riter != str.rend(); riter++) {
        new_str += *riter;
    }
    return new_str;
}

std::string get_complement(std::string& pattern) {
    std::string complement;
    for (auto rev_iter = pattern.rbegin(); rev_iter != pattern.rend();
        rev_iter++) {
        char temp = *rev_iter;
        switch (temp) {
        case 'A':
            complement += 'T';
            break;
        case 'C':
            complement += 'G';
            break;
        case 'G':
            complement += 'C';
            break;
        case 'T':
            complement += 'A';
            break;
        }
    }
    return complement;
}

int get_hamming_distance(std::string str1, std::string str2) {
    int hd = 0;
    for (int i = 0; i < str1.size(); i++) {
        if (str1[i] != str2[i]) hd++;
    }
    return hd;
}

std::vector<std::string> perm_with_repetition(std::vector<std::string> vec,
    std::string nuc, int num) {
    if (num == 0) {
        return vec;
    }
    else if (vec.size() == 0) {
        return perm_with_repetition({ "A", "C", "G", "T" }, nuc, num - 1);
    }
    else {
        std::vector<std::string> vec_new;
        for (auto elem : vec) {
            for (auto& chr : nuc) {
                std::string temp = elem + chr;
                vec_new.push_back(temp);
            }
        }
        vec = std::move(vec_new);
        return perm_with_repetition(vec, nuc, num - 1);
    }
}


// BA01H: find all approximate occurrences of a pattern
std::vector<int> find_all_approx_occur(std::string pattern, std::string genome,
    int num_d) {
    std::vector<int> indices;
    for (int i = 0; i < genome.size() - pattern.size() + 1; i++) {
        std::string frag = genome.substr(i, pattern.size());
        if (get_hamming_distance(pattern, frag) <= num_d) {
            indices.push_back(i);
        }
    }
    return indices;
}

// BA01J: most frequent words with mismatches reverse complement
std::vector<std::string> most_freq_words_with_mismatch_rc(std::string genome,
    int num_k,
    int num_d) {
    int max_freq = 0;
    std::unordered_map<std::string, int> dict;
    std::vector<std::string> kmers = perm_with_repetition({}, "ATGC", num_k);
    for (auto kmer : kmers) {
        int freq1 = find_all_approx_occur(kmer, genome, num_d).size();
        int freq2 =
            find_all_approx_occur(get_complement(kmer), genome, num_d).size();
        int freq = freq1 + freq2;
        if (freq > max_freq) max_freq = freq;
        dict[kmer] = freq;
    }

    std::vector<std::string> answer;
    for (auto iter = dict.begin(); iter != dict.end(); iter++) {
        if (iter->second == max_freq) {
            answer.push_back(iter->first);
        }
    }
    return answer;
}