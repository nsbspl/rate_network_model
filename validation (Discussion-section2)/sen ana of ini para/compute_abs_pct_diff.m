% compute the mean absolute pct diff between two vecs

ref_para = opt_para;
new_para = para_vec_after_44_iter;
pct_diff_vec = abs(new_para - ref_para)./ref_para;
mean_abs_pct_diff_final = mean(pct_diff_vec);