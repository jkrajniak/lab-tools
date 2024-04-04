fmount - easy way to work with sshfs
======================

Requirements
-----------------
- sshfs

Usage
----------
1. copy `fmount` to your `bin/` path
2. put settings in `~/.sshfs`
3. To mount remote directory into local path, type `fmount <local path>`

Settings
-----------

The list of remote paths should be placed in `~/.sshfs` file. The format is quite simple.
You can place default options in `DEFAULT` variable at the beginning of the file. Next, in the following lines you can define the remote paths and the location in local file system.

```bash
ssh hostname:remote_path options local_path
```

For example, the following lines will allow accessing remote `/home/user/backup` directory located on the server `host.pl` by local path `/home/alek/backup`. We set options to `ro` which means that the access will be read-only.

```bash
user@host.pl:/home/user/backup ro /home/alek/backup
```
