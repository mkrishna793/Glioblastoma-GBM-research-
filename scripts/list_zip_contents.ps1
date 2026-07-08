Get-ChildItem -Path 'D:\Glioblastoma Multiforme GBM new research files' -Filter *.zip | ForEach-Object {
    Write-Output "========================================"
    Write-Output "Zip File: $($_.Name)"
    Write-Output "========================================"
    Add-Type -AssemblyName System.IO.Compression.FileSystem
    $zip = [System.IO.Compression.ZipFile]::OpenRead($_.FullName)
    $zip.Entries | Select-Object -First 10 -Property Name, Length
    $zip.Dispose()
}
