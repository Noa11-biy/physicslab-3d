import os, subprocess, sys, tempfile, shutil
SOFFICE = r"C:\Program Files\LibreOffice\program\soffice.exe"
src = os.path.abspath(sys.argv[1])
outdir = os.path.dirname(src)
profile = os.path.join(tempfile.gettempdir(), "lo_profile_cours")
cmd = [SOFFICE, f"-env:UserInstallation=file:///{profile.replace(os.sep, '/')}", "--headless", "--convert-to", "pdf", "--outdir", outdir, src]
r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
print(r.stdout.strip(), r.stderr.strip()[:500])
