# cse284_project
cse284 Project, Spring 2021

## Set up  
From your command line in the folder you want to put the repository:  
`git clone https://github.com/sarah-n-wright/cse284_project.git`

If you don't want to enter your password every time see: https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh to set up ssh key.  

## quick git commands workflow  
`git pull` # pull the latest updates from github  
`git branch <branch name>` # create a new branch  
`git checkout <branch name>` # go to branch <branch name>
`git merge origin/master` # to update your branch to look like main (LB)
 
#### Make your changes on your branch, then:  
`git add <filename>`  # add file to list of files to commit changes  
`git commit -m <description>` # commit all changes to files (that have been been added)  
#### Optional:  
`git push origin <mybranch>`  # make your branch available without merging to main
#### Go back to main `git checkout main` , then: 
`git pull` # pull the latest updates from github in case there have been updates while you worked    
`git merge main <mybranch>` # Merge changes from main to <mybranch>. If this produces conflicts you will get an error message which will need to be fixed.
 
`git merge <mybranch> main` # Merge changes from <mybranch> to main
  
`git push origin main` # push the latest updates from your local main branch to github    
  
## Other useful things:
`git status` # check the status of repository  
`git log` # print a log of recent commits  
`git branch` # show all available branches  

## Pushing other branches  

It can also be useful to push your personal branches:  
`git push origin <mybranch>` as this will make your changes available to everyone, even if they are not ready for merging into main.  
To pull a different branch use:  
`git pull <otherbranch>`
