Adding to this Site
=====================

Basics
-------------------------
Everything here is generated using Sphinx, which is an auto-documentation
program orginally devloped for the Python docs. To build updated versions of
these documentation you will need to install it.

`Install Sphinx <https://www.sphinx-doc.org/en/1.6/install.html>`_

However you don't have to, to start editing. To jump in just navigate to the
sphinx directory in the repo and open up one of the .rst files. .rst files are
what sphinx builds html from which are ultimately rendered to make the website
good pretty.

However you don't have to, to start editing. To jump in just navigate to the
sphinx directory in the repo and open up one of the .rst files. .rst files are
what sphinx builds html from which are ultimately rendered to make the website
good pretty. 

Rst has some special conventions for doing fancy things like adding links or codeblocks.
You can use `this cheatsheet link <https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_
for that kind of stuff. (I have basically just been copy and pasting from here)

Then once you have made some edits, commit them to the repo; very cool!

Rendering rst
-------------
If you want to render the rst to html you will also need to install the Read
the Docs theme. You can do so with the code below.
.. code-block:: bash

   pip install sphinx_rtd_theme

After you have made some changes to rst files navigate to the sphinx directory
in the repo and run :code:`make html`. This will make html files in the
docs/html directory of the repo.

However, becuase I am no expect there are still some jank elements to actually
getting the html to render on the website. You will need to copy all the html
files so they are in the docs dir not docs/html (GitHub likes it better this
way I guess). Then delete the html directory and push your changes to the
repo. To do this run these commands from the CoRNonCOB directory.
.. code-block:: bash

   cp -r ./docs/html/ ./docs/.
   rm -r ./docs/html


Testing
-------
If you want to test your changes before pushing, just do the same make process
and is described above and click on the html file. It should open up in your
browser.



