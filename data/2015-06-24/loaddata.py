import requests
from urllib.error import HTTPError, URLError

def i_filename(i=0, element="CaIR"):
    if element == "Ha" and i<= 811 and i>= 300:
        return f'crisp_l2_20150624_134006_6563_r00{i}.fits'
    elif element == 'CaIR':
        return f'crisp_l2_20150624_134006_8542_r00{i}.fits'
    else:
        raise(f'Given wrong parameters {i=}, {element=}')

def download_files(url):
  """
  Downloads all files from the specified URL.

  Args:
      url: The base URL of the directory containing the files.

  Returns:
      None
  """
  try:
    # Send GET request to retrieve directory listing (potentially)
    response = requests.get(url, allow_redirects=True)
    response.raise_for_status()  # Raise exception for non-2xx status codes

    # # Check for HTML content (may indicate successful directory listing)
    # if response.headers['Content-Type'].startswith('text/html'):
    #   # Extract file names from HTML (implementation may vary depending on website structure)
    #   # This part requires parsing the HTML content to find file links
    #   # For simplicity, this example omits the HTML parsing logic.
    #   # Consider using libraries like Beautiful Soup for robust HTML parsing.
    #   print("Website seems to require HTML parsing for file listing. Implement logic to extract file names from HTML content.")
    #   return

    # Iterate over retrieved content (assuming directory listing format)
    for i in range(589, 812,1):
      # Extract potential filename (adapt based on directory listing format)
      filename = i_filename(i) # Assuming filename is the last element

      # Construct download URL
      download_url = f"{url}/{filename}"

      # Download the file
      download_file(download_url)

  except (HTTPError, URLError) as error:
    print(f"Error downloading files: {error}")

def download_file(url):
  """
  Downloads a single file from the specified URL.

  Args:
      url: The URL of the file to download.

  Returns:
      None
  """
  try:
    response = requests.get(url, stream=True)
    response.raise_for_status()  # Raise exception for non-2xx status codes

    # Get filename from URL (consider potential path handling)
    filename = url.split('/')[-1]

    # Open file for writing in binary mode
    with open(filename, 'wb') as f:
      for chunk in response.iter_content(1024):
        f.write(chunk)

    print(f"Downloaded: {filename}")

  except (HTTPError, URLError) as error:
    print(f"Error downloading file {url}: {error}")

# Replace with the actual URL
base_url = "https://star.pst.qub.ac.uk/webdav/public/fchroma/2015-06-24/8542/"

# Download files from the base URL
download_files(base_url)

print("Download complete (if no errors encountered).")
