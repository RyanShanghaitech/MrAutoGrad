cp -r ../MrAutoGrad_Private/mrautograd_src ./
cp -r ../MrAutoGrad_Private/example ./
cp ../MrAutoGrad_Private/LICENSE ./
cp ../MrAutoGrad_Private/MANIFEST.in ./
cp ../MrAutoGrad_Private/pyproject.toml ./
cp ../MrAutoGrad_Private/README.md ./
cp ../MrAutoGrad_Private/setup.py ./
cp ../MrAutoGrad_Private/install.bash ./
cp ../MrAutoGrad_Private/.gitignore ./

rm -r ./mrautograd_src/ext/mtg # remove the baseline method to respect the copyright