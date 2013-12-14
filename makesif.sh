#!/bin/bash
PASSWORD="password-"$(date +%Y%m%d)
KEY=$1

echo "[Data]" > winnt.sif
echo "AutoPartition = 1" >> winnt.sif
echo "MsDosInitiated = 0" >> winnt.sif
echo "UnattendedInstall = Yes" >> winnt.sif
echo >> winnt.sif
echo "[Unattended]" >> winnt.sif
echo "UnattendMode = FullUnattended" >> winnt.sif
echo "OemSkipEula = Yes" >> winnt.sif
echo "OemPreinstall = No" >> winnt.sif
echo "TargetPath = \WINDOWS" >> winnt.sif
echo "Repartition = Yes" >> winnt.sif
echo "UnattendSwitch = Yes" >> winnt.sif
echo "DriverSigningPolicy = Ignore" >> winnt.sif
echo "WaitForReboot = No" >> winnt.sif
echo >> winnt.sif
echo "[GuiUnattended]" >> winnt.sif
echo "AdminPassword = \"$PASSWORD\"" >> winnt.sif
echo "AutoLogon = Yes" >> winnt.sif
echo "AutoLogonCount = 1" >> winnt.sif
echo "OEMSkipRegional = 1" >> winnt.sif
echo "OemSkipWelcome = 1" >> winnt.sif
echo "TimeZone = 85" >> winnt.sif
echo "ServerWelcome = No" >> winnt.sif
echo >> winnt.sif
echo "[UserData]" >> winnt.sif
echo "ProductKey = \"$KEY\"" >> winnt.sif
echo "FullName = \"user\"" >> winnt.sif
echo "OrgName = \"\"" >> winnt.sif
echo "ComputerName = \"computer\"" >> winnt.sif
echo >> winnt.sif
echo "[Identification]" >> winnt.sif
echo "JoinWorkgroup = WORKGROUP" >> winnt.sif
echo >> winnt.sif
echo "[Networking]" >> winnt.sif
echo "InstallDefaultComponents = Yes" >> winnt.sif
echo >> winnt.sif
echo "[Display]" >> winnt.sif
echo "XResolution = 800" >> winnt.sif
echo "YResolution = 600" >> winnt.sif
echo "BitsPerPel = 32" >> winnt.sif
echo >> winnt.sif
echo "[GuiRunOnce]" >> winnt.sif
echo "Command0=\"E:\VBoxWindowsAdditions.exe /S /with_autologon\"" >> winnt.sif
echo "Command1=\"shutdown -f -s -t 00\"" >> winnt.sif
