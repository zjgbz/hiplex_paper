
conda activate next_cutntag

python ../hiplex_paper/gene_expression_prediction_model/gene_expression_prediction_model_find_params.py

# Define random seeds array
random_seeds=(0 1 7 42 123 999 1234 1337 2021 31415)

# parralel
for current_seed in "${random_seeds[@]}"; do
    echo "Running job with random_seed=$current_seed"
    
    python ../hiplex_paper/gene_expression_prediction_model/gene_expression_prediction_model_train_and_pdp.py \
        $current_seed &
done

wait
echo "All jobs finished!"

