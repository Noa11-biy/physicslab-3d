# Lance PhysicsLab 3D, clique éventuellement à des positions données, prend une capture d'écran en pleine
# résolution et ferme l'application. Sert à vérifier visuellement une simulation sans intervention manuelle.
#
#   .\tools\screenshot.ps1 -Level 5 -Sim 4 -Wait 8
#   .\tools\screenshot.ps1 -Level 5 -Sim 2 -Clicks "176,522" -Wait 3     # coche une case avant la capture
#
# La capture est écrite dans $OutDir (dossier temporaire par défaut), jamais dans le dépôt.
# Les coordonnées de -Clicks sont en pixels physiques de l'écran ("x,y;x,y").
param(
    [int]$Level = 3,
    [int]$Sim = 1,
    [string]$Clicks = "",
    [int]$Wait = 5,
    [string]$Out = "shot.png",
    [string]$OutDir = (Join-Path $env:TEMP "physicslab_shots")
)

$repo = Split-Path -Parent $PSScriptRoot
$exe = Join-Path $repo "build\physicslab.exe"
if (-not (Test-Path $exe)) { throw "Exécutable introuvable : $exe (compiler d'abord : cmake --build build)" }
New-Item -ItemType Directory -Force $OutDir | Out-Null
Remove-Item (Join-Path $OutDir "imgui.ini") -ErrorAction SilentlyContinue  # disposition par défaut à chaque lancement

Add-Type -AssemblyName System.Windows.Forms, System.Drawing
Add-Type @"
using System.Runtime.InteropServices;
public class PhysicsLabNative {
  [DllImport("user32.dll")] public static extern bool SetProcessDPIAware();
  [DllImport("user32.dll")] public static extern bool SetCursorPos(int x, int y);
  [DllImport("user32.dll")] public static extern void mouse_event(uint f, uint dx, uint dy, uint d, int e);
}
"@
# Sans cela, PowerShell ne voit qu'une partie de l'écran sur un écran mis à l'échelle (125 %) : capture tronquée.
[PhysicsLabNative]::SetProcessDPIAware() | Out-Null

# Dossier de travail = $OutDir : l'imgui.ini généré n'atterrit pas dans le dépôt.
$p = Start-Process -FilePath $exe -ArgumentList "--level", $Level, "--sim", $Sim -WorkingDirectory $OutDir -PassThru
Start-Sleep -Seconds 3

foreach ($c in ($Clicks -split ";")) {
    if ($c -eq "") { continue }
    $xy = $c -split ","
    [PhysicsLabNative]::SetCursorPos([int]$xy[0], [int]$xy[1]) | Out-Null
    Start-Sleep -Milliseconds 300
    [PhysicsLabNative]::mouse_event(2, 0, 0, 0, 0)
    Start-Sleep -Milliseconds 120
    [PhysicsLabNative]::mouse_event(4, 0, 0, 0, 0)
    Start-Sleep -Milliseconds 700
}
Start-Sleep -Seconds $Wait

$bounds = [System.Windows.Forms.Screen]::PrimaryScreen.Bounds
$bmp = New-Object System.Drawing.Bitmap $bounds.Width, $bounds.Height
$gfx = [System.Drawing.Graphics]::FromImage($bmp)
$gfx.CopyFromScreen($bounds.Location, [System.Drawing.Point]::Empty, $bounds.Size)
$path = Join-Path $OutDir $Out
$bmp.Save($path, [System.Drawing.Imaging.ImageFormat]::Png)
$gfx.Dispose(); $bmp.Dispose()
Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue
Write-Output $path
