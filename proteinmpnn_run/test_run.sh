conda activate /home/tiejun/miniconda3/envs/mpnn_cu12.4 

path_to_PDB="../pro/6m0j.pdb"
path_to_fasta="../data/6M0J_AE.fasta"
path_to_output="./results"

if [ ! -d $path_to_output ]
then
    mkdir -p $path_to_output
fi

chains_to_design="A E"

python ../../ProteinMPNN/protein_mpnn_run.py \
        --pdb_path $path_to_PDB \
        --pdb_path_chains "$chains_to_design" \
        --out_folder $path_to_output \
        --num_seq_per_target 5 \
        --sampling_temp "0.1" \
        --score_only 1 \
        --seed 0 \
        --batch_size 1
