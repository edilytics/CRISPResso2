# c2eb

### Instructions for installing an instance of CRISPRessoWeb on AWS Elastic Beanstalk.



## Prerequesites

In order to deploy C2Web on AWS Elastic Beanstalk with this guide you will need the following:

- An AWS account with access to EC2, EFS, VPC, ECR, and Elastic Beanstalk
- [AWS EB CLI](https://docs.aws.amazon.com/elasticbeanstalk/latest/dg/eb-cli3.html)
- [Docker](https://docs.docker.com/get-docker/)
- The `c2web` and `c2web-apache` docker images

## AWS EB INSTRUCTIONS

1. Create a Key pair

   - [New Key pair](https://us-east-2.console.aws.amazon.com/ec2/v2/home?region=us-east-2#CreateKeyPair)
   - Set a Name (to be used in step 4)
   - save the .ppk file (for use with PuTTY on Windows) or .pem (for use with OpenSSH on Linux or OSX)

2. Create Security group

   - [New Security Group](https://us-east-2.console.aws.amazon.com/vpc/home?region=us-east-2#SecurityGroups:sort=groupId)
   - Create one for EFS access
   - Set outbound rules: 
     - Type: “All TCP” Destination: “Anywhere”
   - Set inbound rules:
     - Type “NFS” Source: “Anywhere”

3. Create EFS

   - [Create a new filesystem](https://us-east-2.console.aws.amazon.com/efs/home?region=us-east-2#/filesystems)
   - 30 day lifecycle policy - all others default
   - ‘Manage Network Access’ -> add security group in step 1 to all three mount targets
   - [storage-efs-mountfilesystem.config](.ebextensions/storage-efs-mountfilesystem.config) <- add id of filesystem in 2 here

4. Create EB

   - [New Environment](https://us-east-2.console.aws.amazon.com/elasticbeanstalk/home?region=us-east-2#/newEnvironment)
   - **Environment Tier**:Select `Web server environment`
   - **Application information**: Set desired Application name
   - **Environment information**: Set desired Environment name
   - **Platform**: Select `Managed Platform` > Platform: `Docker`, Platform branch: `Running on 64bit Amazon Linux` (not Linux 2 <- important in 2020 when Linux 2 had a bug for mounting file systems)
   - **Application Code**: Select `Sample application`
   - **Service access**: EC2 key pair: AWS security key from step 1
   - **Instances**: security (Add security group from step 1)

5. Upload C2Web image to your ECR private repo

   - [ECR](https://us-east-2.console.aws.amazon.com/ecr/private-registry/repositories)
   - click *Create repository*
   - **General settings**: Select *Private*
   - Repository name: `c2web` (or whatever you prefer)
   - click *Create repository*
   - [Back to ECR](https://us-east-2.console.aws.amazon.com/ecr/private-registry/repositories)
   - click your new repository name
   - click *View push commands* and follow instructions to build images locally and push to your ECR repo

7. Repeat step 5. for the `c2web-apache` docker image

8. Update docker image URIs

   - Update the `docker-compose.yml` file to point to your ECR `c2web` image:
   ```
   web:
      image: &web-image-uri
            '<your-registry-uri>' # replace this with your ECR image URI for the c2web image
   ```
      - Update the `docker-compose.yml` file to point to your ECR `c2web-apache`:
   ```
   apache:
      image: <apache-image-uri-on-ecr> # replace this with your ECR image URI for the c2web-apache image
   ```

7. Deploy (via [EB CLI](https://docs.aws.amazon.com/elasticbeanstalk/latest/dg/eb-cli3.html))

   - `cd` into the project directory
   - run `eb init` -> Select region and application
   - run `eb deploy`

8. Set up forwarding

   - Add zones on bluehost (don’t add subdomain)
     - [https://my.bluehost.com/cgi/dm/zoneedit](https://my.bluehost.com/cgi/dm/zoneedit)
     - Add zone for {demo} -> {elasticbeanstalk.com} (type = CNAME)
     - Add zone for www.demo -> elasticbeanstalk.com (type=CNAME)
   - Add DNS on AWS
     - [https://us-east-2.console.aws.amazon.com/acm/home?region=us-east-2#/wizard/](https://us-east-2.console.aws.amazon.com/acm/home?region=us-east-2#/wizard/)
     - Request a Certificate
     - Type in ‘{demo}.edilytics.com’
     - ‘Add another name’
     - ‘www.{demo}.edilytics.com’
     - Use DNS verification
   - Back in bluehost, add DNS verification CNAMEs
     - (enter once without underscores, then go down and modify and enter underscores)
     - (remove trailing periods)

9. Customize Environment Variables

   - Set environment variables for the C2Web application in the `common-variables` anchor located in `docker-compose.yml`:
   ```
   x-common-variables: &common-variables
   DEBUG: 0
   USE_DEFAULT_USER: 'False'
   ALLOW_ADMIN_INIT: 'True'
   MAMBA_DOCKERFILE_ACTIVATE: 1
   MAX_CONTENT_LENGTH: 2147483648
   FLASK_APP: 'CRISPRessoWEB'
   BANNER_TEXT: 'Look at my cool banner'
   ```

10. Set up SSL (Optional)

   - [https://us-east-2.console.aws.amazon.com/acm/home?region=us-east-2#/certificates/request](https://us-east-2.console.aws.amazon.com/acm/home?region=us-east-2#/certificates/request)
   - Select `Request a public certificate`
   - **Domain Names** Enter the domain name(s) from which you will access your instace (e.g. "crispresso.my-biotech-company.com")
   - **Validation method** select DNS validation if possible
   - Click *Request*
   - Click on the certificate you just requested
   - In your DNS Settings, add the CNAME record shown under **Domains**
   - [Return to EB](https://us-east-2.console.aws.amazon.com/elasticbeanstalk/home?region=us-east-2#/environments)
   - select the environment you created in step 4
   [//]: # (Warning: Instructions below seem to be out-of-date)
   - **Configuration** (in left side bar under the name of environment)
   - **Load balancer** 'Edit' (if there is no 'Edit' button on 'Load Balancer', select **Capacity** 'Edit' >> Environment Type: choose 'Load balanced')
   - **Listeners** 'Add listener'
      - click *Add listener*
      - set Listener port: 443
      - Listener protocol: HTTPS
      - Instance port: 80
      - Instance protocol: HTTP
      - select the certificate then'Save' and 'Apply'