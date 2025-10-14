#!/bin/bash

SFTP_HOST="gweusftp.azenta.com"
SFTP_USER=""
SFTP_PASSWORD=""
FILELIST="files_need_download.txt"
REMOTE_DIR="/90-889954109/00_fastq"

# Use a loop to read each line from the file and download the corresponding file
while IFS= read -r FILENAME
do
    sftp -oBatchMode=no -b - "$SFTP_USER@$SFTP_HOST" <<EOF
        cd $REMOTE_DIR
        get $FILENAME
EOF
done < $FILELIST
