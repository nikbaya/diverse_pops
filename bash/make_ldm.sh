n_cpu=$( grep -c ^processor /proc/cpuinfo )

n_snps=$( cat not_EUR.chr22.maf_gt_0.bim | wc -l )

max_chunk_idx=$((${n_cpu}))

chunk_size=$((${n_snps}/${max_chunk_idx}))


function max16 {
   while [ `jobs | wc -l` -ge 16 ]
   do
      sleep 5
   done
}

function make_ldm {
	local chunk_idx=$1
	
	min_snp_idx=$(( ${chunk_idx}*${chunk_size} + 1 ))
	
	if [ ${chunk_idx} -eq ${max_chunk_idx} ]; then
		max_snp_idx=${n_snps}
	else
		max_snp_idx=$(( (${chunk_idx}+1)*${chunk_size}))
	fi
	
	
	echo "${min_snp_idx}-${max_snp_idx}"
	ldm_type="sparse"

	./gctb \
		--bfile not_EUR.chr22.maf_gt_0 \
		--snp ${min_snp_idx}-${max_snp_idx} \
		--make-${ldm_type}-ldm \
		--out not_EUR.chr22.maf_gt_0
}





for idx in `seq 0 $((${max_chunk_idx}-1))`; do
   max16;  make_ldm ${idx} & 
done
