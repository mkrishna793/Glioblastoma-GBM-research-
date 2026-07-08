import os
import pandas as pd
from alphagenome.models import dna_client

# Set API Key directly in environment for the test
os.environ["ALPHA_GENOME_API_KEY"] = "AIzaSyAtkWXKygAyFOJpCV8CVkXqL8b0z25TvAI"
os.environ["ALPHAGENOME_API_KEY"] = "AIzaSyAtkWXKygAyFOJpCV8CVkXqL8b0z25TvAI"

print("Initializing AlphaGenome client...")
try:
    api_key = os.environ.get("ALPHA_GENOME_API_KEY")
    # Try creating client with default address or dns address as per SKILL.md
    dna_model = dna_client.create(api_key=api_key, address='dns:///gdmscience.googleapis.com:443')
    print("Client created successfully!")
    print("Client properties:", [p for p in dir(dna_model) if not p.startswith('_')])
except Exception as e:
    print(f"Error initializing client: {e}")
