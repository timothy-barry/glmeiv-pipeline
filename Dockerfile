# base image: samsondsc/odm_image
FROM samsondsc/odm_image:1.0

# Install glmeiv from Github
RUN Rscript -e 'devtools::install_github("timothy-barry/glmeiv@HEAD")'
