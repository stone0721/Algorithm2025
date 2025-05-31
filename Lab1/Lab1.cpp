#include <iostream>
#include <vector>
#include <string>
#include <map>
using namespace std;

// DNA序列处理相关函数
string getComplementSeq(const string& seq) {
    map<char, char> pairMap = {
        {'A', 'T'}, {'T', 'A'}, 
        {'C', 'G'}, {'G', 'C'}
    };
    string compSeq;
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        compSeq += pairMap[*it];
    }
    return compSeq;
}

// 计算两个序列的最长公共前缀长度
vector<vector<int>> calcCommonPrefix(const string& seq1, const string& seq2) {
    int len1 = seq1.length(), len2 = seq2.length();
    vector<vector<int>> lenTable(len1 + 1, vector<int>(len2 + 1, 0));
    
    for (int i = len1 - 1; i >= 0; --i) {
        for (int j = len2 - 1; j >= 0; --j) {
            if (seq1[i] == seq2[j]) {
                lenTable[i][j] = lenTable[i+1][j+1] + 1;
            }
        }
    }
    return lenTable;
}

// 计算序列内部的最长公共前缀
vector<vector<int>> calcSelfPrefix(const string& seq, int maxLen) {
    int seqLen = seq.length();
    vector<vector<int>> lenTable(seqLen + 1, vector<int>(seqLen + 1, 0));
    
    for (int i = seqLen - 1; i >= 0; --i) {
        int endPos = (i + maxLen > seqLen - 1) ? seqLen - 1 : i + maxLen;
        for (int j = endPos; j > i; --j) {
            if (seq[i] == seq[j]) {
                lenTable[i][j] = lenTable[i + 1][j + 1] + 1;
            }
        }
    }
    return lenTable;
}

// 查找重复片段
struct RepeatInfo {
    int startPos;
    int endPos;
    string segment;
    int repeatCount;
};

vector<RepeatInfo> findRepeats(const string& seq, int maxLen,
    const vector<vector<int>>& forwardTable,
    const vector<vector<int>>& reverseTable) {
    
    int seqLen = seq.length();
    auto selfTable = calcSelfPrefix(seq, maxLen);
    vector<RepeatInfo> results;
    
    for (int pos = 0; pos < seqLen; ) {
        int bestLen = 0;
        RepeatInfo bestRepeat;
        bool found = false;
        
        // 找出最大匹配长度
        int maxMatchLen = 0;
        for (int j = 0; j < maxLen + 1; ++j) {
            int currLen = forwardTable[pos][j] > reverseTable[pos][j] ? forwardTable[pos][j] : reverseTable[pos][j];
            if (currLen > maxMatchLen) maxMatchLen = currLen;
        }
        
        // 检查不同长度的重复
        for (int len = 1; len <= maxMatchLen; ++len) {
            string unit = seq.substr(pos, len);
            int count = 1;
            
            // 统计重复次数
            while (pos + count * len + len <= seqLen && 
                   selfTable[pos + (count - 1) * len][pos + count * len] >= len) {
                count++;
            }
            
            if (count >= 2) {
                int totalLen = count * len;
                if (totalLen > bestLen) {
                    bestLen = totalLen;
                    bestRepeat = {pos, pos + totalLen, unit, count};
                    found = true;
                }
            }
        }
        
        if (found) {
            results.push_back(bestRepeat);
            pos = bestRepeat.endPos;
        } else {
            pos++;
        }
    }
    
    return results;
}

int main() {
    // 输入序列
    string refSeq = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA";
    string querySeq = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA";
    
    int refLen = refSeq.length();
    string compRefSeq = getComplementSeq(refSeq);
    
    // 计算匹配表
    auto forwardTable = calcCommonPrefix(querySeq, refSeq);
    auto reverseTable = calcCommonPrefix(querySeq, compRefSeq);
    
    // 查找重复
    auto repeats = findRepeats(querySeq, refLen, forwardTable, reverseTable);
    
    // 输出结果
    cout << "重复片段结果:\n";
    cout << "序号\t位置\t长度\t重复次数\t方向\n";
    
    int refIndex = 0, queryIndex = 0, resultIndex = 0;
    
    for (const auto& repeat : repeats) {
        if (repeat.endPos > queryIndex + forwardTable[queryIndex][refIndex]) {
            int tempEnd = repeat.endPos;
            int stopPos = queryIndex + forwardTable[queryIndex][refIndex];
            int repeatCount = 0;
            
            while (tempEnd > stopPos && tempEnd > repeat.startPos) {
                repeatCount++;
                tempEnd -= repeat.segment.length();
            }
            
            // 确定匹配位置和类型
            int searchRange = refSeq.length() - repeat.segment.length();
            int matchPos = -1;
            string matchType;
            
            for (int i = 0; i < searchRange; ++i) {
                if (forwardTable[repeat.startPos][i] >= repeat.segment.length()) {
                    matchType = "正向";
                    matchPos = i + repeat.segment.length();
                    break;
                } 
                if (reverseTable[repeat.startPos][i] >= repeat.segment.length()) {
                    matchType = "反向";
                    matchPos = refSeq.length() - i;
                    break;
                }
            }
            
            resultIndex++;
            cout << resultIndex << "\t" 
                << matchPos << "\t" 
                << repeat.segment.length() << "\t" 
                << repeatCount << "\t\t" 
                << matchType << "\n";
            if (stopPos > repeat.startPos) {
                refIndex += (tempEnd - queryIndex);
            }
            queryIndex = repeat.endPos;
        }
    }
    
    return 0;
}