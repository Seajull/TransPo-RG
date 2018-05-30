from subprocess import check_output
import sys, os

def updateTag() :
    ver=".version"
    if not os.path.isdir(".git"):
        sys.exit("This does not appear to be a Git repository.")
    if not os.path.isfile(".version") :
        out = check_output(["git", "tag"]).decode("utf-8")
        vers=""
        with open(".version","a") as verOut :
            print("File '.version' created.")
            for i in out :
                if i != "v" and i != "\n" :
                    vers+=i
                if i == "\n" :
                    verOut.write(vers+"\n")
                    vers=""
    else :
        out = check_output(["git", "describe","--long", "--tags"]).decode("utf-8")
        vers=(out.split("-")[0][1:]+"\n")
        with open(".version","a+") as verOut :
            verOut.seek(0)
            f=verOut.readlines()
            if f[-1] != vers :
                verOut.write(vers)
            else :
                print("Version file already up-to-date.")

def getVersion() :
    with open(".version","r") as verF :
        f=verF.readlines()
        print("\n"+"Version : "+f[-1])
    return

if __name__ == "__main__" :
    if sys.argv[1] == "update" :
        updateTag()
    if sys.argv[1] == "version" :
        getVersion()
