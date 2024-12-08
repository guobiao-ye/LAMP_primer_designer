import pandas as pd
import joblib
import ast


# 1. Load the trained random forest model and make predictions
def load_and_predict(input_csv, model_path):
    # Load the data
    data = pd.read_csv(input_csv)
    X = data.iloc[:, 1:]  # Features start from the second column

    # Load the trained model
    try:
        model = joblib.load(model_path)
        print("Model loaded successfully.")
    except FileNotFoundError:
        raise FileNotFoundError(f"Model file not found at {model_path}. Please check the path.")
    except Exception as e:
        raise RuntimeError(f"An error occurred while loading the model: {e}")

    # Predict probabilities
    try:
        probabilities = model.predict_proba(X)[:, 1]  # Extract probability of the positive class
    except AttributeError:
        raise AttributeError("The loaded model does not support 'predict_proba'. Ensure it is a scikit-learn model.")

    data['score'] = probabilities
    sorted_data = data.sort_values(by='score', ascending=False)  # Sort by probability
    return sorted_data[['lamp_id', 'score']]


# 2. Extract and unpack detailed information from full_scores
def extract_full_scores(input_csv):
    data = pd.read_csv(input_csv)
    extracted_data = []

    for _, row in data.iterrows():
        lamp_id = row['lamp_id']
        try:
            full_scores = ast.literal_eval(row['full_scores'])  # Convert string to dictionary
            single_primer_scores = full_scores['single_primer_scores']
            distances = full_scores['distances']
            dimer_score = full_scores['dimer_score']
        except (ValueError, SyntaxError):
            raise ValueError(f"Invalid 'full_scores' format for lamp_id: {lamp_id}")

        # Initialize the result row
        result_row = {'lamp_id': lamp_id, 'dimer_score': dimer_score}

        # Extract information for each primer
        for primer in single_primer_scores:
            primer_name = primer['primer']
            result_row[f"{primer_name}_sequence"] = primer['sequence']
            result_row[f"{primer_name}_gc"] = primer['gc']
            result_row[f"{primer_name}_tm"] = primer['tm']
            result_row[f"{primer_name}_delta_g_5p"] = primer['delta_g_5p']
            result_row[f"{primer_name}_delta_g_3p"] = primer['delta_g_3p']
            result_row[f"{primer_name}_hairpin"] = primer['hairpin']
            result_row[f"{primer_name}_gc_score"] = primer['gc_score']
            result_row[f"{primer_name}_tm_score"] = primer['tm_score']
            result_row[f"{primer_name}_delta_g_score"] = primer['delta_g_score']
            result_row[f"{primer_name}_hairpin_score"] = primer['hairpin_score']

        # Extract distance information
        for distance_name, value in distances.items():
            result_row[distance_name] = value

        extracted_data.append(result_row)

    return pd.DataFrame(extracted_data)


# 3. Merge probabilities and detailed information
def merge_results(probabilities_csv="data/output/Intermediate_file/score_into_model.csv",
                  full_scores_csv="data/output/Intermediate_file/full_score_results.csv",
                  model_path="model/random_forest_model.pkl",
                  output_csv="data/output/Final_score/final_sorted_results.csv"):
    try:
        probabilities = load_and_predict(probabilities_csv, model_path)
        full_scores = extract_full_scores(full_scores_csv)
    except Exception as e:
        raise RuntimeError(f"An error occurred during the merging process: {e}")

    # Merge the two parts of data
    merged_data = pd.merge(probabilities, full_scores, on='lamp_id')
    merged_data.to_csv(output_csv, index=False)


# Main program
if __name__ == "__main__":
    # Input file paths
    probabilities_input = "data/output/Intermediate_file/score_into_model.csv"
    full_scores_input = "data/output/Intermediate_file/full_score_results.csv"
    model_path = "model/random_forest_model.pkl"
    output_file = "data/output/Final_score/final_sorted_results.csv"

    # Run the merging process
    try:
        merge_results(probabilities_input, full_scores_input, model_path, output_file)
        print(f"The final results have been saved to {output_file}")
    except Exception as e:
        print(f"An error occurred: {e}")
