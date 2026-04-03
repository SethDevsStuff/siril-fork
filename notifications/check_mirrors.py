import json
import requests
import random
import os

def check_mirrors(json_file):
    # 1. Load the JSON data
    if not os.path.exists(json_file):
        print(f"Error: {json_file} not found.")
        return

    with open(json_file, 'r') as f:
        mirrors = json.load(f)

    # Configuration
    file_prefix = "siril_cat1_healpix8_xpsamp_"
    file_suffix = ".dat"
    num_files = 48  # 0 through 47
    checks_per_file = 1
    chunk_size = 100 # bytes

    for mirror in mirrors:
        base_url = mirror['url'].rstrip('/')
        description = mirror['description']
        
        print(f"\n{'='*60}")
        print(f"TESTING MIRROR: {description}")
        print(f"BASE URL: {base_url}")
        print(f"{'='*60}")

        for i in range(num_files):
            filename = f"{file_prefix}{i}{file_suffix}"
            
            # Note: Some mirrors like Zenodo might have different path structures
            # but we follow your format: base_url + filename
            full_url = f"{base_url}/{filename}"
            
            print(f"Checking {filename}...")
            
            for _ in range(checks_per_file):
                # Using a safe range (0-5000) for initial bytes
                start_byte = random.randint(0, 5000)
                end_byte = start_byte + chunk_size
                
                headers = {
                        'Range': f'bytes={start_byte}-{end_byte}',
                        'User-Agent': 'Siril/0.0-mirror-check'
                        }
                
                try:
                    # stream=True is more efficient for just checking the response code/headers
                    response = requests.get(full_url, headers=headers, timeout=10, stream=True)
                    
                    if response.status_code == 206:
                        print(f"  [PASS] Range {start_byte}-{end_byte}: 206 Partial Content")
                    else:
                        print(f"  [WARN] {filename}: Status {response.status_code}")
                        sys.exit(1)
                    
                    # Close connection early since we only need the header/status
                    response.close()
                        
                except requests.exceptions.RequestException as e:
                    print(f"  [ERR] Connection failed: {e}")
                    break # Skip to next file if connection is dead

if __name__ == "__main__":
    check_mirrors('spcc_mirrors.json')
