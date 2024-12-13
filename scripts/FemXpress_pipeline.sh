python FemXpress.py \
	-b '/data5/wangxin/20220923fgh_mouse_brain/mouse_brain_cellranger/220731C_mousebrain/outs1/possorted_genome_bam.bam' \
	-g '/data2/wangxin/database/genome/mm10.fa' \
	-e '/data5/wangxin/20220923fgh_mouse_brain/mouse_brain_cellranger/meta/mb_brain_meta.txt' \
	-r '/data2/wangxin/database/processed_genome/mm10.rmsk.txt' \
	-n 'mouse' \
	-a '/data2/wangxin/database/genome/gencode.vM25.chr_patch_hapl_scaff.basic.annotation.gtf' \
	-c 2
