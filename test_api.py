#!/usr/bin/env python3
"""
Test script to check ESM Atlas API connectivity and endpoints.
Run this to diagnose API issues before using the main app.
"""

import requests
import time
import sys

def test_endpoint(url, sequence="MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"):
    """Test a single API endpoint."""
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    
    try:
        print(f"Testing: {url}")
        response = requests.post(url, headers=headers, data=sequence, timeout=30)
        print(f"  Status: {response.status_code}")
        
        if response.status_code == 200:
            print(f"  ✅ Success! Response length: {len(response.content)} bytes")
            return True
        elif response.status_code == 503:
            print(f"  ⚠️  Service temporarily unavailable")
        elif response.status_code == 429:
            print(f"  ⚠️  Rate limited")
        else:
            print(f"  ❌ Error: {response.text[:100]}...")
            
    except requests.exceptions.Timeout:
        print(f"  ⏱️  Timeout")
    except requests.exceptions.ConnectionError:
        print(f"  🔌 Connection error")
    except Exception as e:
        print(f"  ❌ Error: {e}")
    
    return False

def main():
    """Test all available ESM Atlas API endpoints."""
    print("🧪 Testing ESM Atlas API Endpoints")
    print("=" * 50)
    
    endpoints = [
        'https://api.esmatlas.com/foldSequence/v1/pdb/',
        'https://esmatlas.com/foldSequence/v1/pdb/',
        'https://api.esmatlas.com/foldSequence/v1/pdb'
    ]
    
    working_endpoints = []
    
    for endpoint in endpoints:
        if test_endpoint(endpoint):
            working_endpoints.append(endpoint)
        print()
        time.sleep(1)  # Be nice to the API
    
    print("=" * 50)
    if working_endpoints:
        print(f"✅ Found {len(working_endpoints)} working endpoint(s):")
        for endpoint in working_endpoints:
            print(f"  - {endpoint}")
    else:
        print("❌ No working endpoints found. The API may be down.")
        print("\nTroubleshooting tips:")
        print("1. Check your internet connection")
        print("2. Wait a few minutes and try again")
        print("3. The ESM Atlas API may be under maintenance")
        print("4. Consider using local ESMFold instead")
    
    return len(working_endpoints) > 0

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
