from openai import OpenAI
import openai
import os
import json
import time
import datetime

#%%
##### INITIAL TARGET GUESSES FROM NCI #####
# Initialize client
openai.api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

# Paths
input_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_data/"
output_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_targets/"

# Model to use
MODEL = "gpt-4o"

# Wait time after each chunk (in seconds)
WAIT_AFTER_CHUNK = 90

# Process all chunks
for i in range(1,220): 
    filename = f"nci_drug_definitions_chunk_{i}_non_batch.jsonl"
    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, f"response_chunk_{i}.jsonl")

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
                record_id = record.get("id")

                # Submit to chat completion endpoint
                response = client.chat.completions.create(
                    model=MODEL,
                    messages=messages,
                    temperature=0,
                )

                responses.append({
                    "id": record_id,
                    "response": response.choices[0].message.content
                })

                time.sleep(1.1)  # avoid rate limits

            except Exception as e:
                print(f"Error processing a record: {e}")
                responses.append({
                    "id": record.get("id"),
                    "response": None,
                    "error": str(e)
                })

    # Save to output JSONL
    with open(output_path, "w") as out_f:
        for r in responses:
            out_f.write(json.dumps(r) + "\n")

    print(f"Finished chunk {i} of 219 and saved results to {output_path}")
    time.sleep(WAIT_AFTER_CHUNK)

#%%
##### REFINE TARGET GUESSES FROM NCI #####
openai.api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

# Paths
input_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_data/"
output_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_targets/"

# Model to use
MODEL = "gpt-4o"

# Wait time after each chunk (in seconds)
WAIT_AFTER_CHUNK = 90 

# Process all chunks
for i in range(1,220):
    filename = f"nci_drug_definitions_refine_chunk_{i}_non_batch.jsonl"
    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, f"response_refine_chunk_{i}.jsonl")

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
                record_id = record.get("id")

                # Submit to chat completion endpoint
                response = client.chat.completions.create(
                    model=MODEL,
                    messages=messages,
                    temperature=0,
                )

                responses.append({
                    "id": record_id,
                    "response_refine": response.choices[0].message.content
                })

                time.sleep(1.1)  

            except Exception as e:
                print(f"Error processing a record: {e}")
                responses.append({
                    "id": record.get("id"),
                    "response_refine": None,
                    "error": str(e)
                })

    # Save to output JSONL
    with open(output_path, "w") as out_f:
        for r in responses:
            out_f.write(json.dumps(r) + "\n")

    print(f"Finished chunk {i} of 219 and saved results to {output_path}")
    time.sleep(WAIT_AFTER_CHUNK)
 
#%%
##### INITIAL TARGET GUESSES FROM DRUGS #####
openai.api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

# Paths
input_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_data/"
output_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_targets/"

# Model to use
MODEL = "gpt-4o"

# Wait time after each chunk (in seconds)
WAIT_AFTER_CHUNK = 90 

# Process all chunks
for i in range(1,26): 
    filename = f"clinical_trial_drugs_chunk_{i}_non_batch.jsonl"
    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, f"response_drugs_chunk_{i}.jsonl")

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
                record_id = record.get("drug")

                # Submit to chat completion endpoint
                response = client.chat.completions.create(
                    model=MODEL,
                    messages=messages,
                    temperature=0,
                )

                responses.append({
                    "drug": record_id,
                    "response_drugname": response.choices[0].message.content
                })

                time.sleep(1.1)  

            except Exception as e:
                print(f"Error processing a record: {e}")
                responses.append({
                    "drug": record.get("drug"),
                    "response_drugname": None,
                    "error": str(e)
                })

    # Save to output JSONL
    with open(output_path, "w") as out_f:
        for r in responses:
            out_f.write(json.dumps(r) + "\n")

    print(f"Finished chunk {i} of 25 and saved results to {output_path}")
    time.sleep(WAIT_AFTER_CHUNK)
    
#%%
##### REFINE TARGET GUESSES FROM DRUGS #####
openai.api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

# Paths
input_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_data/"
output_dir = "/mnt/Data/Callahan/projects/trial_selector/data/openai_chunks_targets/"

# Model to use
MODEL = "gpt-4o"

# Wait time after each chunk (in seconds)
WAIT_AFTER_CHUNK = 90 

# Process all chunks
for i in range(1,26):
    filename = f"clinical_trial_drugs_refine_chunk_{i}_non_batch.jsonl"
    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, f"response_drugs_refine_chunk_{i}.jsonl")

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
                record_id = record.get("drug")

                # Submit to chat completion endpoint
                response = client.chat.completions.create(
                    model=MODEL,
                    messages=messages,
                    temperature=0,
                )

                responses.append({
                    "drug": record_id,
                    "response_refine": response.choices[0].message.content
                })

                time.sleep(1.1)  # avoid rate limits

            except Exception as e:
                print(f"Error processing a record: {e}")
                responses.append({
                    "drug": record.get("drug"),
                    "response_refine": None,
                    "error": str(e)
                })

    # Save to output JSONL
    with open(output_path, "w") as out_f:
        for r in responses:
            out_f.write(json.dumps(r) + "\n")

    print(f"Finished chunk {i} of 25 and saved results to {output_path}")
    time.sleep(WAIT_AFTER_CHUNK)
