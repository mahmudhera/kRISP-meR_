from on_target_scoring.rs2_score_calculator import get_rs2_score
import sys
import os

model_file = 'V3_model_nopos.pickle'
score = get_rs2_score('G'*30, model_file)
print(score)
score = get_rs2_score('G'*30, model_file)
print(score)
