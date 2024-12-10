* Add the files `Fonts` and `geoviewer.so` to the path `flair-2.3`;
* Running, using: 
```
python flair.py
```






### PUSH a new branch to a exiting `git`:

#### Connecting a Local File to the Repository
1. **Initializing the Local Repository (if not already initialized)**
   - If your local file is not yet associated with any Git repository, first open the terminal (Command Prompt or PowerShell on Windows, Terminal application on Linux or macOS) in the directory where the local file is located. Then, run the `git init` command. This will create a hidden `.git` folder in the local directory to manage the version control information of the local repository.

2. **Adding the Remote Repository Association**
   - Open the terminal and navigate to the directory of the local file. Use the command `git remote add origin https://github.com/xinwenir/flair.git`. Here, `origin` is an alias for the remote repository. You can use other names, but `origin` is commonly used. This command establishes an association between the local repository and the remote repository `https://github.com/xinwenir/flair.git`.

3. **Adding Local Files to the Local Repository and Committing**
   - Use the command `git add .` to add all files in the local directory to Git's staging area. If you only want to add specific files, you can use `git add [filename]` (for example, `git add myfile.txt`).
   - Run the command `git commit -m "Initial commit"` (where `"Initial commit"` is the commit message and can be modified according to the actual situation). This commits the files in the staging area to the local repository, recording the version changes of the files locally.
4. **Pushing Local Files to the Remote Repository**
   - If the remote repository has protected branches (such as the `master` or `main` branch), you may need to first allow your account to push in the settings of the remote repository or create a new branch for pushing.
   - Run the command `git push -u origin [branch name]`. If your local repository currently has only one branch (usually `master` or `main`), you can directly use `git push -u origin master` (if it is the `main` branch, use `git push -u origin main`). The `-u` option sets the upstream branch, so that you can simply use the `git push` command to push updates later.

```
git init
git remote add origin https://github.com/xinwenir/flair.git
git checkout -b flair-dev-PY
git add . # (or: git add [filename])
git commit -m "Initial commit"
git push -u origin flair-dev-PY
```