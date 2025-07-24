from openai import OpenAI
import openai
import os
import json
import time

#%%
##### Summary #####
# Initialize client
openai.api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

# Paths
input_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_therapies_curated_negative/"
output_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_therapies_curated_negative_answers/"

# Model to use
MODEL = "gpt-4o"

# Wait time after each chunk (in seconds)
WAIT_AFTER_CHUNK = 90 # so we don't hit TPM limit

# Process all chunks
for i in range(1,22): # 21 chunks per guess method
    filename = f"clinical_trial_therapy_guess_summary_chunk_{i}.jsonl"
    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, f"response_clinical_trial_therapy_guess_summary_chunk_{i}.jsonl")

    if os.path.exists(output_path):
        print(f"Skipping {filename}, already processed.")
        continue

    print(f"Processing {filename}...")
    responses = []

    with open(input_path, "r") as f:
        for line in f:
            try:
                record = json.loads(line)
                messages = record["messages"]
                record_id = record.get("trial")

                # Submit to chat completion endpoint
                response = client.chat.completions.create(
                    model=MODEL,
                    messages=messages,
                    temperature=0,
                )

                responses.append({
                    "trial": record_id,
                    "response": response.choices[0].message.content
                })

                time.sleep(1.1)  # avoid rate limits

            except Exception as e:
                print(f"Error processing a record: {e}")
                responses.append({
                    "trial": record.get("trial"),
                    "response": None,
                    "error": str(e)
                })

    # Save to output JSONL
    with open(output_path, "w") as out_f:
        for r in responses:
            out_f.write(json.dumps(r) + "\n")

    print(f"Finished chunk {i} of 21 and saved results to {output_path}")
    time.sleep(WAIT_AFTER_CHUNK)

#%%
##### Title #####
# Initialize client
openai.api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

# Paths
input_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_therapies_curated_negative/"
output_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_therapies_curated_negative_answers/"

# Model to use
MODEL = "gpt-4o"

# Wait time after each chunk (in seconds)
WAIT_AFTER_CHUNK = 90 # so we don't hit TPM limit

# Process all chunks
for i in range(1,22): # 21 chunks per guess method
    filename = f"clinical_trial_therapy_guess_title_chunk_{i}.jsonl"
    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, f"response_clinical_trial_therapy_guess_title_chunk_{i}.jsonl")

    if os.path.exists(output_path):
        print(f"Skipping {filename}, already processed.")
        continue

    print(f"Processing {filename}...")
    responses = []

    with open(input_path, "r") as f:
        for line in f:
            try:
                record = json.loads(line)
                messages = record["messages"]
                record_id = record.get("trial")

                # Submit to chat completion endpoint
                response = client.chat.completions.create(
                    model=MODEL,
                    messages=messages,
                    temperature=0,
                )

                responses.append({
                    "trial": record_id,
                    "response": response.choices[0].message.content
                })

                time.sleep(1.1)  # avoid rate limits

            except Exception as e:
                print(f"Error processing a record: {e}")
                responses.append({
                    "trial": record.get("trial"),
                    "response": None,
                    "error": str(e)
                })

    # Save to output JSONL
    with open(output_path, "w") as out_f:
        for r in responses:
            out_f.write(json.dumps(r) + "\n")

    print(f"Finished chunk {i} of 21 and saved results to {output_path}")
    time.sleep(WAIT_AFTER_CHUNK)

 
#%%
##### Intervention #####
# Initialize client
openai.api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

# Paths
input_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_therapies_curated_negative/"
output_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_therapies_curated_negative_answers/"

# Model to use
MODEL = "gpt-4o"

# Wait time after each chunk (in seconds)
WAIT_AFTER_CHUNK = 90 # so we don't hit TPM limit

# Process all chunks
for i in range(1,22): # 21 chunks per guess method
    filename = f"clinical_trial_therapy_guess_intervention_chunk_{i}.jsonl"
    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, f"response_clinical_trial_therapy_guess_intervention_chunk_{i}.jsonl")

    if os.path.exists(output_path):
        print(f"Skipping {filename}, already processed.")
        continue

    print(f"Processing {filename}...")
    responses = []

    with open(input_path, "r") as f:
        for line in f:
            try:
                record = json.loads(line)
                messages = record["messages"]
                record_id = record.get("trial")

                # Submit to chat completion endpoint
                response = client.chat.completions.create(
                    model=MODEL,
                    messages=messages,
                    temperature=0,
                )

                responses.append({
                    "trial": record_id,
                    "response": response.choices[0].message.content
                })

                time.sleep(1.1)  # avoid rate limits

            except Exception as e:
                print(f"Error processing a record: {e}")
                responses.append({
                    "trial": record.get("trial"),
                    "response": None,
                    "error": str(e)
                })

    # Save to output JSONL
    with open(output_path, "w") as out_f:
        for r in responses:
            out_f.write(json.dumps(r) + "\n")

    print(f"Finished chunk {i} of 21 and saved results to {output_path}")
    time.sleep(WAIT_AFTER_CHUNK)
    