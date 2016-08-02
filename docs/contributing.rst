.. currentmodule:: nlabpy
.. _contributing:

************
Contributing
************

Contributions are welcome!


Pull Requests
~~~~~~~~~~~~~

Some guidelines on contributing to ggplot:

* Please submit a Pull Requests.
* Pull Requests can be submitted as soon as there is code worth discussing.
  Pull Requests track the branch, so you can continue to work after the PR is submitted.
  Review and discussion can begin well before the work is complete,
  and the more discussion the better.
  The worst case is that the PR is closed.
* Pull Requests should generally be made against master
* Pull Requests should be tested, if feasible:
    - bugfixes should include regression tests
    - new behavior should at least get minimal exercise
    - see tests in `nlabpy/tests/` subdir and the information above.


Some git
~~~~~~~~

This assumes that you cloned from the original ggplot repository (`origin`) and added your
own fork (`mine`) like this:

.. code-block:: bash

    git clone https://github.com/eco32i/nlabpy.git # becomes "origin"
    git remote add mine git@github.com:YourName/ggplot.git


* Use a feature branch: `git checkout origin/master ; git checkout -b <feature_name>`
* Commit edits belonging together in one go: `git add <files>; git commit`
* If you need to add or change something to your last commit: `git commit --amend`
* Push to your own clone: `git push mine` (during the first time it will tell you to
  use some more arguments to setup the push for this branch)
* Add commits to your PR: just push them : `git commit; git push mine`
* Rebase against master get the latest changes in master and to see if your changes
  still apply cleanly : `git rebase origin/master`. If there are some conflicting
  changes you will be prompted to cleanup the merge conflict.
* Do not use `git merge origin/master` or your PR will be cluttered with merge 
  messages. (Note: this is only valid for your PRs, not if you push your own code
  to your own repo for others to work against!)
* Force push after a rebase: `git push -f mine`.
* Make changes to old commits (reorder, edits, squash multiple commits into 
  one): `git rebase -i origin/master`:
  
  - reorder the commits be moving lines up and down
  - prefix the commit you want to edit with 'e' and the ones you want to squash
    into the former one with 's'
  - for the edits: use `git commit --amend` or `git commit` and after all your
    edits are done: `git rebase --continue`
  - for the squash: just edit the commit message

* View the commit history: `git log`. To only view the first line of the commit
  message: `git log --oneline`
* view the diff of already `git add`ed files: `git diff --cached`
* Unadd files: `git reset HEAD -- <filename>`. You can add an 'unadd' alias for this
  via `git config --global alias.unadd "reset HEAD --"` -> `git unadd <file>`
