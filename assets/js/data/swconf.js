<<<<<<< HEAD
<<<<<<< HEAD
const swconf = {
  
    cacheName: 'chirpy-1733275908',resources: [
      '/assets/css/jekyll-theme-chirpy.css',
      '/',
      
        '/categories/',
      
        '/tags/',
      
        '/archives/',
      
        '/about/',
      

      
      
    ],

    interceptor: {paths: [
        
      ],urlPrefixes: [
        
=======
---
layout: compress
permalink: '/:path/swconf.js'
# Note that this file will be fetched by the ServiceWorker, so it will not be cached.
---

=======
>>>>>>> main
const swconf = {
  
    cacheName: 'chirpy-1733363963',resources: [
      '/assets/css/jekyll-theme-chirpy.css',
      '/',
      
        '/categories/',
      
        '/tags/',
      
        '/archives/',
      
        '/about/',
      

      
      
    ],

<<<<<<< HEAD
    interceptor: {
      {%- comment -%} URLs containing the following paths will not be cached. {%- endcomment -%}
      paths: [
        {% for path in site.pwa.cache.deny_paths %}
          {% unless path == empty %}
            '{{ path | relative_url }}'{%- unless forloop.last -%},{%- endunless -%}
          {% endunless  %}
        {% endfor %}
      ],

      {%- comment -%} URLs containing the following prefixes will not be cached. {%- endcomment -%}
      urlPrefixes: [
        {% if site.analytics.goatcounter.id != nil and site.pageviews.provider == 'goatcounter' %}
          'https://{{ site.analytics.goatcounter.id }}.goatcounter.com/counter/'
        {% endif %}
>>>>>>> main
=======
    interceptor: {paths: [
        
      ],urlPrefixes: [
        
>>>>>>> main
      ]
    },

    purge: false
<<<<<<< HEAD
<<<<<<< HEAD
  
};

=======
  {% else %}
    purge: true
  {% endif %}
};
>>>>>>> main
=======
  
};

>>>>>>> main
