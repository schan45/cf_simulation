$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot
$python = if ($env:TDK_PYTHON) { $env:TDK_PYTHON } else { 'python' }
& $python -m streamlit run app.py --server.address localhost --browser.gatherUsageStats false
exit $LASTEXITCODE
