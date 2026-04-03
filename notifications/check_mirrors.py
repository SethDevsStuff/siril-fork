import json
import requests
import random
import sys
import os

def check_mirrors(json_file):
    if not os.path.exists(json_file):
        print(f"Error: {json_file} not found.")
        sys.exit(1)

    try:
        with open(json_file, 'r') as f:
            mirrors = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON: {e}")
        sys.exit(1)

    # Configuration
    file_prefix = "siril_cat1_healpix8_xpsamp_"
    file_suffix = ".dat"
    num_files = 48
    checks_per_file = 1
    chunk_size = 100 

    # Using a Session object keeps the connection alive
    session = requests.Session()
    session.headers.update({'User-Agent': 'Siril/1.5-mirror-check'})

    for mirror in mirrors:
        base_url = mirror['url'].rstrip('/')
        description = mirror['description']
        
        print(f"\n--- TESTING MIRROR: {description} ({base_url}) ---")

        for i in range(num_files):
            filename = f"{file_prefix}{i}{file_suffix}"
            full_url = f"{base_url}/{filename}"
            
            for _ in range(checks_per_file):
                start_byte = random.randint(0, 5000)
                end_byte = start_byte + chunk_size
                
                # Range header is specific to each request
                headers = {'Range': f'bytes={start_byte}-{end_byte}'}
                
                try:
                    # stream=False (default) + Session for connection reuse
                    response = session.get(full_url, headers=headers, timeout=10)
                    
                    # Read the response body into memory
                    data = response.content
                    
                    if response.status_code == 206:
                        # Verify we actually got the number of bytes we expected
                        if len(data) == (chunk_size + 1):
                            print(f"  [PASS] Range {start_byte}-{end_byte}: 206 Partial Content")
                            continue # Pass
                        else:
                            print(f"[FAILURE] {filename}: Byte count mismatch. Got {len(data)} bytes.")
                            sys.exit(1)
                    else:
                        print(f"[FAILURE] {full_url}")
                        print(f"Expected 206, but received {response.status_code} for file {filename}")
                        #sys.exit(1)
                        
                except requests.exceptions.RequestException as e:
                    print(f"[CONNECTION FAILURE] {full_url}")
                    print(f"Error details: {e}")
                    sys.exit(1)

        print(f"Mirror '{description}' passed all checks.")

    print("\nAll mirrors verified successfully.")
    sys.exit(0)

if __name__ == "__main__":
    check_mirrors('spcc_mirrors.json')
