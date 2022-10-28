
import hashlib
def get_file_md5(filename):
    file_hash = hashlib.md5()
    with open(filename, "rb") as f:
        while chunk := f.read(8192):
            file_hash.update(chunk)
    return  file_hash.hexdigest()
