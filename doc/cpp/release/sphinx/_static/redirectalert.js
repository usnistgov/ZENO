function redirectalert(url) {
	if (confirm("Thank you for visiting NIST.\n\nWe have provided a link to this site because it has information that may be of interest to our users. NIST does not necessarily endorse the views expressed or the facts presented on this site. Further, NIST does not endorse any commercial products that may be advertised or available on this site.\n\nClick OK to be directed to: "+url)) {
        	window.location = url
	}
}
