#%%
from openai import OpenAI
import openai
import os
import json
import time

#%%
##### CURATED TARGET CELL SURFACE #####
# Initialize client
openai.api_key = os.getenv("OPENAI_API_KEY")
client = OpenAI()

# Paths
input_dir = "/mnt/Data/Callahan/projects/trial_selector/data/raw/"
output_dir = "/mnt/Data/Callahan/projects/trial_selector/data/processed/"

# Model to use
MODEL = "gpt-4o"

# Wait time after each chunk (in seconds)
WAIT_AFTER_CHUNK = 90 # so we don't hit TPM limit

# loop overjsonl 3 times
for i in range(1,4): 
    filename = "cell_surface_curated.jsonl"
    input_path = os.path.join(input_dir, filename)
    output_path = os.path.join(output_dir, f"cell_surface_curated_response_{i}.jsonl")

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

    print(f"Finished replicate {i} of 3 and saved results to {output_path}")
    time.sleep(WAIT_AFTER_CHUNK)

