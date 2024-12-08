import pandas as pd  
import joblib  

# 1. Load the model  
model_path = 'model/random_forest_model.pkl'  
rf_model = joblib.load(model_path)  

# 2. Load the new dataset for prediction  
# Replace 'path_to_your_prediction_data.csv' with the actual path to your CSV file  
data_to_predict = pd.read_csv('data/primers_scoring_output.csv')  

# 3. Prepare the features for prediction  
# Extract features (excluding the 'lamp_id' column)  
X_new = data_to_predict.iloc[:, 1:]  

# 4. Predict probabilities for the normal primers (label=1)  
predicted_probabilities = rf_model.predict_proba(X_new)[:, 1]  # Get the probabilities of label 1  

# 5. Create a DataFrame to show lamp_ids and their corresponding probabilities  
results = pd.DataFrame({  
    'lamp_id': data_to_predict['lamp_id'],  
    'probability_of_normal_primers': predicted_probabilities  
})  

# 6. Sort the results by probability in descending order  
sorted_results = results.sort_values(by='probability_of_normal_primers', ascending=False).reset_index(drop=True)  

# 7. Display sorted results  
print(sorted_results)  

# Optional: Save sorted results to a CSV file  
sorted_results.to_csv('predicted_probabilities.csv', index=False)