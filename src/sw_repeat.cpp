#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
using namespace Rcpp;

inline int idx(int i, int j, int ncols) {
  return i * ncols + j;
}



// [[Rcpp::export(rng=FALSE)]]
List sw_circular_cpp2(std::string sequence1,std::string sequence2, bool output_matrix=false,
                     int match_score = 1,
                     int mismatch_penalty = -1,
                     int gap_penalty = -2,
                     int diagonal_penalty = -10000,
                     int offset = 0) {
  
  std::string seq_str1 = (sequence1);
  std::string seq_str2 = (sequence2);
  if(seq_str1.length()!=seq_str2.length() ){
    return( List::create() );
  }
  
  int n = seq_str1.length();
  //offset = n / 2-std::min(n/4, 250);
  //offset = n / 2;
  
  // 计算每行需要存储的列数
  // 注意：原代码只填充 j < i-offset 的部分
  int cols_needed=1;
  std::vector<int> row_sizes(n + 1, 0);
  for (int i = offset; i <= n; i++) {
    //int cols_needed = std::max(0, i - offset);
    row_sizes[i] = cols_needed++;
    
  }
  
  // 使用压缩存储：每行只存需要的列（j从1开始）
  std::vector<std::vector<int>> score_mat(n + 1);
  std::vector<std::vector<char>> trace_mat(n + 1);
  
  for (int i = offset; i <= n; i++) {
    score_mat[i].resize(row_sizes[i] + 1, 0);  // +1 因为j从1开始，索引0不用
    trace_mat[i].resize(row_sizes[i] + 1, 's');
  }
  
  // 设置对角线惩罚（对应原代码中预先设置的部分）
  // 原代码：for (int i = offset; i <= n; i++) { score_mat(i, jj) = diagonal_penalty; jj++; }
  // 这里设置的是 score_mat[i][i-offset] 位置（因为 jj 从0开始，对应 j = i-offset）
  int jj = 0;
  
  for (int i = offset; i <= n; i++) {
    // if (j < row_sizes[i]) {  // 确保索引有效
    //   score_mat[i][j] = diagonal_penalty;
    //   // 注意：这里不设置trace，因为原代码也没设置trace
    // }
    score_mat[i][jj]=diagonal_penalty;
    jj++;
    
  }
  
  int start_max_i=0;
  int start_max_j=0;
  int total_max_score=0;
  
  // 填充下三角矩阵（只填充 j < i-offset 的部分）
  for (int i = offset+1; i <= n; i++) {
    int max_j = i - offset;
    for (int j = 1; j < max_j; j++) {
      // 获取对角得分 (i-1, j-1)
      int diag_score = 0;
      if (j-1 < row_sizes[i-1] && i-1 >= 1) {
        diag_score = score_mat[i-1][j-1];
      }
      
      // 计算匹配/不匹配得分（原代码中已注释掉 i-1 == j-1 的判断）
      if (seq_str1[i-1] == seq_str2[j-1]) {
        diag_score += match_score;
      } else {
        diag_score += mismatch_penalty;
      }
      
      // 上方向得分 (i-1, j)
      int up_score = gap_penalty;
      if (j < row_sizes[i-1] && i-1 >= 1) {
        up_score = score_mat[i-1][j] + gap_penalty;
      }
      
      // 左方向得分 (i, j-1)
      int left_score = gap_penalty;
      if (j-1 < row_sizes[i]) {
        left_score = score_mat[i][j-1] + gap_penalty;
      }
      
      // 取最大值
      int max_score = diag_score;
      char trace = 'd';
      
      if (left_score > max_score) {
        max_score = left_score;
        trace = 'l';
      }
      if (up_score > max_score) {
        max_score = up_score;
        trace = 'u';
      }
      
      score_mat[i][j] = max_score;
      trace_mat[i][j] = trace;
      
      if(max_score>total_max_score){
        total_max_score=max_score;
        start_max_i=i;
        start_max_j=j;
      }
      
      
    }
  }
  
  // 找回溯起点（最后一行最大值）
  int max_score_val = 0;
  int start_j = 1;int start_i=n;
  std::vector<int> last_row_scores;
  
  if (row_sizes[n] > 0) {
    for (int j = 1; j < row_sizes[n]; j++) {
      int score = score_mat[n][j];
      last_row_scores.push_back(score);
      
      if (score > max_score_val) {
        max_score_val = score;
        start_j = j;
      }
    }
  }
  
  
  
  start_j=start_max_j;
  start_i=start_max_i;
  max_score_val=total_max_score;
  
  
  int normal_score=1;
  int normal_score_first50=10000;
  int normal_score_last50=100000;
  
  // 回溯路径
  std::vector<int> path_i, path_j;
  path_i.push_back(start_i);
  path_j.push_back(start_j);
  int i = start_i, j = start_j;
  std::string align1, align2;

  while (i > 0 && j > 0 && j < row_sizes[i]) {
    char trace = trace_mat[i][j];
    if (trace == 's') break;

    path_i.push_back(i-1);
    path_j.push_back(j-1);

    if (trace == 'd') {
      //align1 = seq_str[i-1] + align1;
      //align2 = seq_str[j-1] + align2;
      
      if (seq_str1[i-1] == seq_str2[j-1]) {
        normal_score++;
        
      }else{
        normal_score--;
        
      }  
      
      i--; j--;
    } else if (trace == 'u') {
      //align1 = seq_str[i-1] + align1;
      //align2 = "-" + align2;
      normal_score--;
      
      i--;
    } else if (trace == 'l') {
      //align1 = "-" + align1;
      //align2 = seq_str[j-1] + align2;
      normal_score--;
      
      j--;
    } else {
      break;
    }
    
    if(j==start_j-100-1){
      normal_score_last50=normal_score;
    }
    
    if(j==100+1){
      normal_score_first50=normal_score;
    }
    
  }
  
  normal_score_first50=normal_score-normal_score_first50;
  
  int path_len = (path_i.size()+path_j.size() )/2;
  bool reached_first_col = true;
  double circular_score = (normal_score > 0) ? (double)normal_score / (double)path_len : 0;
  
  
  
  if(!output_matrix){
  
  // Named("is_circular") = reached_first_col &&    
  //   (    (circular_score > 0.25 && path_len >= 200)  || (circular_score > 0.4 && path_len <= 100)  
  //          || (circular_score > 0.3 && path_len < 200 && path_len > 100)   ) &&
  //            start_j >= 15  ,
  
  return List::create(
    Named("is_circular") = reached_first_col &&    
      (    (circular_score >= 0.25 )  ) && ( path_len*1.0/n > 0.125) && 
      start_j >= 15 ,
      Named("score") = circular_score,
      Named("max_score") = max_score_val,
      Named("normal_score") = normal_score,
      Named("normal_score_first50") = normal_score_first50,
      Named("normal_score_last50") = normal_score_last50,
      
      Named("path_length") = path_len,
      Named("reached_first_col") = reached_first_col,
      Named("alignment") = align1 + "\n" + align2,
      Named("alignment1") = align1,
      Named("alignment2") = align2,
      Named("start_j") = start_j,
      Named("start_i") = start_i,
      
      Named("path_i") = wrap(path_i),
      Named("path_j") = wrap(path_j),
      Named("last_row_scores")=wrap(last_row_scores)
    );

  }else{
  
  
   i = start_i, j = start_j;normal_score=1;
    std::string align1, align2;
    
    while (i > 0 && j > 0 && j < row_sizes[i]) {
      char trace = trace_mat[i][j];
      if (trace == 's') break;
      
      path_i.push_back(i-1);
      path_j.push_back(j-1);
      
      if (trace == 'd') {
        align1 = seq_str1[i-1] + align1;
        align2 = seq_str2[j-1] + align2;
        
        if (seq_str1[i-1] == seq_str2[j-1]) {
          normal_score++;
        }else{
          normal_score--;
          
        }  
        
        i--; j--;
      } else if (trace == 'u') {
        //align1 = seq_str[i-1] + align1;
        //align2 = "-" + align2;
        normal_score--;
        
        i--;
      } else if (trace == 'l') {
        align1 = "-" + align1;
        align2 = seq_str2[j-1] + align2;
        normal_score--;
        
        j--;
      } else {
        break;
      }
    }
  

    Rcpp::NumericMatrix score_mat_full(n + 1, n + 1);
    Rcpp::CharacterMatrix trace_mat_full(n + 1, n + 1);
  
    // 初始化trace矩阵为"s"
    for (int r = 0; r <= n; r++) {
      for (int c = 0; c <= n; c++) {
        trace_mat_full(r, c) = "s";
      }
    }
  
    // 将压缩存储的数据复制到完整矩阵中
    for (int r = 1; r <= n; r++) {
      for (int c = 1; c < row_sizes[r]; c++) {
        score_mat_full(r, c) = score_mat[r][c];
        // 转换char到string
        if (trace_mat[r][c] == 'd') {
          trace_mat_full(r, c) = "d";
        } else if (trace_mat[r][c] == 'u') {
          trace_mat_full(r, c) = "u";
        } else if (trace_mat[r][c] == 'l') {
          trace_mat_full(r, c) = "l";
        } else if (trace_mat[r][c] == 's') {
          trace_mat_full(r, c) = "s";
        }
      }
    }
  
  //&& normal_score_last50 >0 && normal_score_first50>0
  
  // Named("is_circular") = reached_first_col &&    
  //   (    (circular_score > 0.25 && path_len >= 200)  || (circular_score > 0.4 && path_len <= 100)  
  //          || (circular_score > 0.3 && path_len < 200 && path_len > 100)   ) &&
  //            start_j >= 15  ,
  
    return List::create(
      Named("is_circular") = reached_first_col &&    
        (    (circular_score >= 0.25 )  ) && ( path_len*1.0/n > 0.125) && 
        start_j >= 15 ,
      Named("score") = circular_score,
      Named("max_score") = max_score_val,
      Named("normal_score") = normal_score,
      Named("normal_score_first50") = normal_score_first50,
      Named("normal_score_last50") = normal_score_last50,
      Named("path_length") = path_len,
      Named("reached_first_col") = reached_first_col,
      Named("alignment") = align1 + "\n" + align2,
      Named("alignment1") = align1,
      Named("alignment2") = align2,
      Named("start_j") = start_j,
      Named("start_i") = start_i,
      
      Named("score_matrix") = score_mat_full,  // 新增：返回完整的得分矩阵
      Named("trace_matrix") = trace_mat_full,  // 新增：返回完整的回溯矩阵
      Named("path_i") = wrap(path_i),
      Named("path_j") = wrap(path_j),
      Named("last_row_scores")=wrap(last_row_scores)
      
    );

  }
  
  return List::create();
  
}



// [[Rcpp::export(rng=FALSE)]]
List sw_circular_cpp(std::string sequence, bool output_matrix=false,
                     int match_score = 1,
                     int mismatch_penalty = -1,
                     int gap_penalty = -3,
                     int diagonal_penalty = -10000,
                     int offset = 0) {
  
  int n=sequence.length();
  int t_offset= floor( n / 2.0-std::min(n/4.0, 2050.0) );
  
  return(sw_circular_cpp2(sequence,sequence,
                          output_matrix=output_matrix,
                          match_score=match_score,
                          mismatch_penalty=mismatch_penalty,
                          gap_penalty=gap_penalty,
                          diagonal_penalty=diagonal_penalty,
                          offset= t_offset
                          
                          ) );
  
  
}







// [[Rcpp::export(rng=FALSE)]]
List sw_cpp(CharacterVector sequence1, 
                      CharacterVector sequence2,
                      int match_score = 1,
                      int mismatch_penalty = -1,
                      int gap_penalty = -3) {
  
  std::string seq1_str = Rcpp::as<std::string>(sequence1);
  std::string seq2_str = Rcpp::as<std::string>(sequence2);
  int n1 = seq1_str.length();
  int n2 = seq2_str.length();
  
  // 预转换字符为索引
  auto char_to_idx = [](char c) {
    switch(c) {
    case 'A': return 0;
    case 'T': return 1;
    case 'C': return 2;
    case 'G': return 3;
    default: return -1;
    }
  };
  
  // 预计算匹配矩阵
  int match_matrix[4][4];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      match_matrix[i][j] = (i == j) ? match_score : mismatch_penalty;
    }
  }
  
  // 转换序列为索引
  std::vector<int> seq1_idx(n1);
  std::vector<int> seq2_idx(n2);
  for (int i = 0; i < n1; i++) seq1_idx[i] = char_to_idx(seq1_str[i]);
  for (int j = 0; j < n2; j++) seq2_idx[j] = char_to_idx(seq2_str[j]);
  
  // 使用一维数组
  int n2p1 = n2 + 1;
  std::vector<int> score_mat((n1 + 1) * n2p1, 0);
  std::vector<unsigned char> trace_mat((n1 + 1) * n2p1, 0);
  
  // 辅助索引函数
  auto idx = [n2p1](int i, int j) { return i * n2p1 + j; };
  
  int max_score = 0;
  int max_i = 0, max_j = 0;
  
  // 主循环
  for (int i = 1; i <= n1; i++) {
    int c1_idx = seq1_idx[i-1];
    int row_offset = i * n2p1;
    int prev_row_offset = (i-1) * n2p1;
    
    for (int j = 1; j <= n2; j++) {
      int c2_idx = seq2_idx[j-1];
      
      int diag_score = score_mat[prev_row_offset + (j-1)] + match_matrix[c1_idx][c2_idx];
      int up_score = score_mat[prev_row_offset + j] + gap_penalty;
      int left_score = score_mat[row_offset + (j-1)] + gap_penalty;
      
      int max_cell = diag_score;
      unsigned char direction = 1;
      
      if (up_score > max_cell) {
        max_cell = up_score;
        direction = 2;
      }
      if (left_score > max_cell) {
        max_cell = left_score;
        direction = 3;
      }
      if (0 > max_cell) {
        max_cell = 0;
        direction = 0;
      }
      
      int current_idx = row_offset + j;
      score_mat[current_idx] = max_cell;
      trace_mat[current_idx] = direction;
      
      if (max_cell > max_score) {
        max_score = max_cell;
        max_i = i;
        max_j = j;
      }
    }
  }
  
  // 优化回溯：使用 vector 追加，最后反转（O(n) 而非 O(n²)）
  std::vector<char> align1_vec, align2_vec;
  align1_vec.reserve(max_score + 10);
  align2_vec.reserve(max_score + 10);
  
  int i = max_i, j = max_j;
  int normal_score = 0;
  
  while (i > 0 && j > 0) {
    unsigned char dir = trace_mat[idx(i, j)];
    if (dir == 0) break;
    
    if (dir == 1) {  // diag
      align1_vec.push_back(seq1_str[i-1]);
      align2_vec.push_back(seq2_str[j-1]);
      if (seq1_str[i-1] == seq2_str[j-1])
        normal_score++;
      else
        normal_score--;
      i--; j--;
    } else if (dir == 2) {  // up
      align1_vec.push_back(seq1_str[i-1]);
      align2_vec.push_back('-');
      normal_score--;
      i--;
    } else if (dir == 3) {  // left
      align1_vec.push_back('-');
      align2_vec.push_back(seq2_str[j-1]);
      normal_score--;
      j--;
    }
  }
  
  // 反转得到正确顺序
  std::reverse(align1_vec.begin(), align1_vec.end());
  std::reverse(align2_vec.begin(), align2_vec.end());
  
  std::string align1(align1_vec.begin(), align1_vec.end());
  std::string align2(align2_vec.begin(), align2_vec.end());
  
  // 提取最后一行的分数
  IntegerVector last_row_score(n2 + 1);
  for (int j = 0; j <= n2; j++) {
    last_row_score[j] = score_mat[idx(n1, j)];
  }
  
  return List::create(
    Named("max_score") = max_score,
    Named("last_row_scores") = last_row_score,
    Named("normal_score") = normal_score,
    Named("path_length") = (int)align1.length(),
    Named("alignment") = align1 + "\n" + align2,
    Named("alignment1") = align1,
    Named("alignment2") = align2,
    Named("max_i") = max_i,
    Named("max_j") = max_j
  );
}


// [[Rcpp::export]]
NumericMatrix compute_matches(CharacterVector seq1, CharacterVector seq2,
                              int window, int step, double min_match) {
  int len1 = seq1.size(), len2 = seq2.size();
  int n1 = (len1 - window) / step + 1;
  int n2 = (len2 - window) / step + 1;
  
  NumericMatrix result(n1, n2);
  
  for(int i = 0; i < n1; i++) {
    for(int j = 0; j < n2; j++) {
      int matches = 0;
      for(int k = 0; k < window; k++) {
        if(seq1[i*step + k] == seq2[j*step + k]) matches++;
      }
      double rate = (double)matches / window;
      if(rate >= min_match) result(i, j) = rate;
      else result(i, j) = 0;
    }
  }
  return result;
}


// [[Rcpp::export(rng=FALSE)]]
bool detect_two_poly_regions_cpp(const std::string& seq_str, const IntegerMatrix& t_pos_matrix, int t_threshold=20) {
  
  if (t_pos_matrix.nrow() == 0) return false;
  
  // 筛选长度>5的区域并合并距离<10的区域
  int curr_start = -1, curr_end = -1;
  int max_len = 0;  // 定义并初始化 max_len
  
  for (int i = 0; i < t_pos_matrix.nrow(); i++) {
    int start = t_pos_matrix(i, 0);
    int end = t_pos_matrix(i, 1);
    int len = end - start + 1;
    
    if (len <= t_threshold/4) continue;  // 跳过短区域（阈值20时跳过长度<=5的区域）
    
    if (curr_start == -1) {
      curr_start = start;
      curr_end = end;
      max_len = std::max(max_len, len);  // 更新最大长度
    } else if (start - curr_end < 5) {
      curr_end = end;  // 合并
      int merged_len = curr_end - curr_start + 1;
      max_len = std::max(max_len, merged_len);  // 更新合并后的长度
    } else {
      if (curr_start > 50 && max_len > t_threshold) {
        //std::cout<< max_len <<" "<< start <<" "<<end<< std::endl;
        // std::cout<< max_len <<" "<< curr_start <<" "<<curr_end<< std::endl;
        
        return true;  // 检查上一个合并区域
        
      }
      
      curr_start = start;
      curr_end = end;
      max_len = len;  // 重置最大长度
    }
  }
  
  //std::cout<< max_len << std::endl;
  // 检查最后一个区域
  return max_len > t_threshold && curr_start > 50;
}
