.. MsSpec-ASE documentation master file, created by
   sphinx-quickstart on Mon Jan 18 10:59:09 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. image:: download_btn.png
    :scale: 0%

.. image:: title.png
    :scale: 0%



.. raw:: html

    <style>
        div#mybanner {
                position:relative;
                width:100%;
                height:400px;
                background:linear-gradient(to top, #004371 0%, #ffffff 100%);}
        div#myfooter {
                width:100%;
                height:20px;
                color:white;
                text-align:right;
                font-size:12px;
                background:black;}
        div#myfooter a:link, div#myfooter a:visited{
                color: white;
                text-decoration: none;
                background-color: transparent;
                }
        div#myfooter a:hover, div#myfooter a:active{
                color: white;
                text-decoration: none;
                background-color: underline;
                }
        div#mybutton {
                width:96px;
                position: absolute;
                right: 32px;
                bottom: 20px}
        div#mytitle {
                position: absolute;
                width: 740px;
                left: 40px;
                bottom: 20px;
                }
        div#mycontacts {
                border-left: thin solid #ffffff8e;
                padding-left: 10px;
                width:200px;
                height:50px;
                color:black;
                text-align:left;
                font-size:14px;
                line-height: 1.5em;
                background:none;
                position: absolute;
                left: 30px;
                bottom: 40px;}
        div#mycontacts a:link, div#mycontacts a:hover, div#mycontacts a:active{
                color: white;
                font-style: italic;
                text-decoration: none;
                }
        div#mycontacts a:hover, div#mycontacts a:active{
                color: white;
                font-weight: bold;
                }
    </style>

    <div id="mybanner">
        <div id="mytitle">
          <img src="./_images/title.png" alt="MsSpec">
        </div>
        <div id="mybutton">
            <a href="javascript:DownloadAndRedirect()">
                <img src="./_images/download_btn.png" alt="Download">
            </a>
        </div>

        <div id="mycontacts">
            <a href="mailto:didier.sebilleau@univ-rennes1.fr">
                Didier Sébilleau
            </a>
            </br>
            <a href="mailto:sylvain.tricot@univ-rennes1.fr">
                Sylvain Tricot
            </a>

        </div>
    </div>




    <div id="myfooter">
        <a href="https://www.ipr.univ-rennes1.fr">
            Rennes Institute of Physics, France
        </a>
    </div>

    <script type="text/javascript">
        function DownloadAndRedirect()
        {
            var DownloadURL = "./downloads.html";
            var RedirectURL = "./downloads.html";
            var RedirectPauseSeconds = 1;
            location.href = DownloadURL;
            setTimeout("DoTheRedirect('"+RedirectURL+"')",parseInt(RedirectPauseSeconds*1000));
        }
        function DoTheRedirect(url) { window.location=url; }
    </script>



#################
Table of Contents
#################

.. toctree::
    :maxdepth: 1

    intro
    downloads
    spectroscopies/spectro
    parameters
    tutorials/index
    contributors
    faq/index
    todo



.. only:: html

    ##################
    Indices and tables
    ##################

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`


.. only:: latex

    #################################
    Application Programming Interface
    #################################

    .. toctree::
        :glob:

        modules/*
