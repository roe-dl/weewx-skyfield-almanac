# Language files

## What are those files for?

This extension provides tags to use in WeeWX skins. Some of those tags
return strings. Some of those strings are language dependent. If you want 
this extension return localized strings, for example for the planets' names, 
the constellations' names, or the Venus and Mercury phases, you need 
appropriate entries in the language files of your skin. The files
in this directory provide the required information. 

## Where to put the data out of the files?

Unfortunately you cannot simply put the files somewhere to get them 
recognized. Instead the language files of your skin will have to be
augmented by the data of the files from this directory. Fortunately
the installation of this extension can do this automatically for skins
that use the standard WeeWX way for internationalization.

Please note that the installation does not overwrite any existing
value in the skin's language files. Missing keys are added only.
Language files, that have HTML code in keys (in keys, not values)
cannot be augmented automatically.

## How to add a localization for your local language?

To add another language get the file `template.conf`, copy it to a
file named after the language code of your choosen language,  and add the
names in that language at the right side of the equal signs.
In case there is no local name for some entry, remove the line
of that entry entirely.

## Links

* [User manual and instructions](../README.md)
* [WeeWX localization instructions](https://weewx.com/docs/latest/custom/localization/)
