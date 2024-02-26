# Computational setup
There are lots of ways to customize your computational environment. From shell aliases that make navigating on the command-line easier, to custom settings for version control, to custom functions in R, here are a few things that I like to set up for every new computer that I work on.

# .bashrc, .bash_profile, .zshrc, .Rprofile, etc.
These are the names of files that exist in your home directory on Unix(-like) systems like Linux and MacOS. The dot at the beginning of the name tells the operating system that they are 'hidden' files - so you won't see them in your file manager (e.g. Finder) or with unix commands like `ls`, unless you use special settings or options.
These files are basically scripts that are automatically run every time you want to use a given code interpreter. This is useful because you can create aliases and variables that act as shortcuts for tasks that you repeat often.
- `.bash_profile` is run every time you use `bash` - whether it's running a script or interactively interacting with your system via the Terminal.
- `.bashrc` is *only* run when you are using `bash` interactively, e.g. when you open Terminal on a Mac.
- `.zshrc` is equivalent to `.bashrc`, but for the `zsh` shell scripting language, instead of the `bash` shell scripting language. This has become the default shell for MacOS.
- `.Rprofile` is equivalent to `.bash_profile`, but for the *R* language. This is run when you open the R app, start a new session in Rstudio, run an R script, or type something like `R` in the Terminal.

Other similar files exist for other languages, or even for different, specific instances in these languages. For example, there are files that are sourced system-wide, those specific to each user, those that are *only* sourced when a script is run, but not for interactive use...

# Work in progress
I'm using the fact that I'm writing this 'guide' as an opportunity to actually standardize my own setup, which has been a bit ad-hoc until now. So, what I *currently* do is not necessarily ideal and is definitely subject to change.

At the moment, I pretty much put all of my shell customizations in `.bash_profile`, then add `source ~/.bash_profile` to my `.bashrc` and `.zshrc`. This defeats the purpose of all the nuances in when they are run, but it makes it simpler for me to remember where my settings are. They're always in `.bash_profile`.

# ssh aliases
I log into the same servers over and over again, all the time. Rather than constantly having to type `ssh username@host.domain` all the time, I add `alias ssh_xx='ssh username@host.domain'` to my `.bash_profile` so all I have to do is type `ssh_xx` as a shortcut.
This is especially useful when the server uses a special port or other options, or when you can use an ssh key instead of a password. The alias uses the appropriate options for the ssh command, which can be lengthy.

# R convenience
I have the line `alias newr='open -n -a R.app'` in my `.bash_profile` so it's easy to open new windows of the R app. I should just use Rstudio more, but when I'm just exploring some quick data informally, I find that too complicated / too much screenspace. I often just open a new terminal window and type 'R' to do things interatively there, instead of the app, though.

I have also added `headc <- function(x,nr=5,nc=nr) { head(x,nr)[,1:nc] }` to my `.Rprofile` so I can easily use `headc` instead of `head` to take a peek at data frames and matrices that might have thousands of columns as well as rows. There are tools to do this within the Tidyverse but for quick exploration I don't like loading all those complicated packages and using special data formats...

I will put more convenient functions for the `.Rprofile` file in `./r_functions/convenience.r`
