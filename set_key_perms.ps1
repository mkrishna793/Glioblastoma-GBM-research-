$path = "D:\research of the GBM\ssh\lightning_rsa"
$acl = Get-Acl $path
$acl.SetAccessRuleProtection($true, $false)
$username = [System.Security.Principal.WindowsIdentity]::GetCurrent().Name
$rule = New-Object System.Security.AccessControl.FileSystemAccessRule($username, 'Read', 'Allow')
$acl.SetAccessRule($rule)
Set-Acl -Path $path -AclObject $acl
Write-Output "Successfully set permissions for $username on $path"
