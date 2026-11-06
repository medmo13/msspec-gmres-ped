##[Ps1 To Exe]
##
##Kd3HDZOFADWE8uK1
##Nc3NCtDXThU=
##Kd3HFJGZHWLWoLaVvnQnhQ==
##LM/RF4eFHHGZ7/K1
##K8rLFtDXTiW5
##OsHQCZGeTiiZ4tI=
##OcrLFtDXTiW5
##LM/BD5WYTiiZ4tI=
##McvWDJ+OTiiZ4tI=
##OMvOC56PFnzN8u+Vs1Q=
##M9jHFoeYB2Hc8u+Vs1Q=
##PdrWFpmIG2HcofKIo2QX
##OMfRFJyLFzWE8uK1
##KsfMAp/KUzWJ0g==
##OsfOAYaPHGbQvbyVvnQX
##LNzNAIWJGmPcoKHc7Do3uAuO
##LNzNAIWJGnvYv7eVvnQX
##M9zLA5mED3nfu77Q7TV64AuzAgg=
##NcDWAYKED3nfu77Q7TV64AuzAgg=
##OMvRB4KDHmHQvbyVvnQX
##P8HPFJGEFzWE8tI=
##KNzDAJWHD2fS8u+Vgw==
##P8HSHYKDCX3N8u+Vgw==
##LNzLEpGeC3fMu77Ro2k3hQ==
##L97HB5mLAnfMu77Ro2k3hQ==
##P8HPCZWEGmaZ7/K1
##L8/UAdDXTlGDjpra7jFL9l/8S2skevm/trWyyYSy6/nQjCzXTZUDWmR4gSzuN0O4Vf4uZvYHvcEFRiEnPOEb57GeHv+sJQ==
##Kc/BRM3KXhU=
##
##
##fd6a9f26a06ea3bc99616d4851b372ba
# PowerShell

$SCRIPT_PATH = $MyInvocation.MyCommand.Path
$IMAGE_NAME  = "iprcnrs/msspec:latest"

Function Read-YesNo($PROMPT, $DEFAULT) {
    $VALID = 1
    While($VALID -ne 0) {
        $ANSWER = Read-Host "$PROMPT (y/n) [$DEFAULT]? "
        Switch($ANSWER) {
            ""      {$ANSWER=$DEFAULT; $VALID=0; Break}
            "y"     {$VALID=0; Break}
            "n"     {$VALID=0; Break}
            Default {Write-Host "Invalid choice, please answer 'y' or 'n'."; $VALID=1; Break}
        }
    }
    Return $ANSWER
}

Function Get-MsSpecContainerID {
    Return docker ps -a -q --filter ancestor=$IMAGE_NAME
}

Function Get-MsSpecImageID {
    Return docker images -q $IMAGE_NAME
}

Function Add-MsSpecImage {
	Write-Host "Pulling MsSpec image...";
	docker pull $IMAGE_NAME;
	Write-Host "done."
	$IID = Get-MsSpecImageID;
	Return $IID;
}

Function Add-MsSpecContainer {
    $CID = Get-MsSpecContainerID
    If($CID -eq $NULL) {
		$IID = Get-MsSpecImageID
		If($IID -eq $NULL) {
			#Add-MsSpecImage;
		}
        Write-Host -nonewline "Creating the MsSpec container... "
        $mydocuments = [environment]::GetFolderPath("mydocuments")
        $drive = ($mydocuments -split ":",2)[0]
        docker create -i -v ${drive}:\:/$drive --entrypoint /bin/bash $IMAGE_NAME >$NULL
        $CID = Get-MsSpecContainerID
		Write-Host "done."
    }
    Return $CID
}

Function Start-MsSpecContainer {
    $CID = Add-MsSpecContainer;
	$status = docker container inspect -f '{{.State.Status}}' $CID;
	IF($status -ne "running") {
		Write-Host -nonewline "Starting MsSpec container...";
		$CID = docker start $CID;
		Write-Host "done.";
	}
}

Function Remove-MsSpecContainer {
	Write-Host -nonewline "Removing MsSpec container...";
    $CID = Get-MsSpecContainerID;
	IF($CID -ne $NULL) {
		docker stop -t 1 $CID >$NULL;
		docker rm $CID >$NULL;
	}
	Write-Host "done."
}

Function Run-MsSpec($opts) {
    # Run msspec command
    $CID = Get-MsSpecContainerID;
    $nixPath = "/" + (${PWD} -replace "\\","/") -replace ":",""
    docker exec -i -t -e DISPLAY=host.docker.internal:0 -w $nixPath $CID msspec $opts
}

Function Run-Bash {
    # Run bash
    $CID = Get-MsSpecContainerID;
    $nixPath = "/" + (${PWD} -replace "\\","/") -replace ":",""
    docker exec -i -t -w $nixPath -e DISPLAY=host.docker.internal:0 $CID /bin/bash
}

Function Uninstall-MsSpecImage {
    $ANSWER = Read-YesNo "You are about to remove the MsSpec Docker image, its container and this script. Are you sure" "n"
    Switch($ANSWER) {
        "y" {Remove-MsSpecContainer;
		     Write-Host -nonewline "Removing $IMAGE_NAME...";
			 $IID = Get-MsSpecImageID;
			 If($IID -ne $NULL) {
				docker rmi $IMAGE_NAME >$NULL;
			 }
			 Write-Host "done.";
#             Remove-Item $THIS_SCRIPT;
             Break;
            }
        "n" {Write-Host "Uninstallation aborted."; Break}
    }
}

Function Show-Help {
    Write-Host "Usage: 1) msspec -p [PYTHON OPTIONS] SCRIPT [ARGUMENTS...]"
    Write-Host "       2) msspec [-l FILE | -i | -h]"
    Write-Host "       3) msspec [bash | reset]"
    Write-Host ""
    Write-Host "Form (1) is used to launch a script"
    Write-Host "Form (2) is used to load a hdf5 data file"
    Write-Host "Form (3) is used to control the Docker container/image."
    Write-Host ""
    Write-Host "List of possible options:"
    Write-Host "    -p   Pass every arguments after this option to the msspec"
    Write-Host "         virtual environment Python interpreter."
    Write-Host "    -i   Run the interactive Python interpreter within msspec"
    Write-Host "         virtual environment."
    Write-Host "    -l   Load and display a *.hdf5 data file in a graphical"
    Write-Host "         window."
    Write-Host "    -v   Print the version."
    Write-Host "    -h   Show this help message."
    Write-Host ""
    Write-Host "    bash        This command starts an interactive bash shell in the"
    Write-Host "                MsSpec container."
    Write-Host "    reset       This command removes the MsSpec container (but not the"
    Write-Host "                image). Changes made in the container will be lost and"
    Write-Host "                any new call to msspec will recreate a new fresh container."
}

Switch($args[0]) {
    "reset" {Remove-MsSpecContainer; Break}
    "bash"  {Start-MsSpecContainer; Run-Bash; Break}
	"-u"    {Uninstall-MsSpecImage; Break}
	"-h"    {Start-MsSpecContainer; Show-Help; Break}
	""      {Start-MsSpecContainer; Show-Help; Break}
    Default {Start-MsSpecContainer; Run-MsSpec $args; Break}
	pull    {Add-MsSpecImage; Break}
}
